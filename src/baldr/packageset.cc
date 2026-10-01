#include "baldr/packageset.h"
#include "baldr/complexrestriction.h"
#include "baldr/graphtile.h"
#include "baldr/tilehierarchy.h"
#include "baldr/tilepatch.h"
#include "config.h"
#include "midgard/constants.h"
#include "midgard/logging.h"
#include "midgard/sequence.h"

#include <boost/property_tree/ptree.hpp>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <condition_variable>
#include <filesystem>
#include <fstream>
#include <limits>
#include <mutex>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <unordered_set>

using namespace valhalla::midgard;

namespace {

// Cell size of the RegionPolygon segment index.
constexpr double kCellDeg = 0.05;
constexpr double kMetersPerDegLat = 110574.0;
constexpr double kMetersPerDegLng = 111320.0;
// Depths closer than this are a tie; the newer build wins.
constexpr double kDepthTieMeters = 0.5;
// Cell size of the fine ownership cache used inside level-2 tiles an ownership boundary crosses.
constexpr double kOwnerCellDeg = 0.01;
// Slack for the local equirectangular projection of distances across a box.
constexpr double kProjectionSlack = 1.05;
// Factor turning a projected distance into a lower bound of the true one.
constexpr double kLowerBoundFactor = 0.95;
// Rings of cells RegionPolygon::Distance searches before it falls back to every segment.
constexpr int32_t kMaxDistanceRings = 200;
// Smallest join tolerance a node grid is built for, in meters.
constexpr double kMinGridCellMeters = 1.0;
// Latitude margin, in degrees, for points probed just outside a tile's node grid.
constexpr double kGridLatMarginDeg = 0.1;
// Padding of the reach box around each node, in degrees (well above float precision).
constexpr float kReachPadDeg = 1e-5f;

// Tile reference index: "VTRF", then its format version.
constexpr uint32_t kTileRefsMagic = 0x46525456;
constexpr uint32_t kTileRefsVersion = 2;
// An overlay is one file per package: "VPKJ", its format version, the join key, the number of
// patched tiles, the number of copy ids, a hash of everything after this 32-byte header, then the
// copy ids the join used, a directory (tile id, offset, size) of the patches in tile id order, and
// the patches themselves (see MakeTilePatch).
constexpr uint32_t kOverlayMagic = 0x4a4b5056;
constexpr uint32_t kOverlayVersion = 2;
constexpr size_t kOverlayTileCountOffset = 16;
constexpr size_t kOverlayIdCountOffset = 20;
constexpr size_t kOverlayHashOffset = 24;
constexpr size_t kOverlayHeaderSize = 32;
constexpr size_t kOverlayIdSize = 16;
constexpr size_t kOverlayDirectorySize = 12;
// zlib level of the patches; higher levels save few bytes
constexpr int kPatchLevel = 9;
// Version of the join itself: a change that changes joined tiles must bump it, so overlays
// written before are stale.
constexpr uint64_t kJoinVersion = 1;
// Bytes of package tiles an ahead-of-time join rewrites between drops of the mapped pages.
constexpr size_t kReleaseBytes = size_t(64) << 20;

// Key of the grid cell of the given size (degrees) holding a coordinate.
uint64_t CellKey(const PointLL& ll, double cell) {
  const auto row = static_cast<int64_t>(std::floor(ll.lat() / cell));
  const auto col = static_cast<int64_t>(std::floor(ll.lng() / cell));
  return (static_cast<uint64_t>(row) << 32) ^ static_cast<uint32_t>(col);
}

// Distance in meters between two nearby points (local equirectangular approximation).
double Meters(const PointLL& a, const PointLL& b) {
  const double kx = kMetersPerDegLng * std::cos((a.lat() + b.lat()) * 0.5 * kPiD / 180.0);
  const double dx = (a.lng() - b.lng()) * kx, dy = (a.lat() - b.lat()) * kMetersPerDegLat;
  return std::sqrt(dx * dx + dy * dy);
}

// Distance in meters from p to segment ab, in a local equirectangular projection around p.
double SegmentDistance(const PointLL& p, const PointLL& a, const PointLL& b) {
  const double kx = kMetersPerDegLng * std::cos(p.lat() * kPiD / 180.0);
  const double ax = (a.lng() - p.lng()) * kx, ay = (a.lat() - p.lat()) * kMetersPerDegLat;
  const double bx = (b.lng() - p.lng()) * kx, by = (b.lat() - p.lat()) * kMetersPerDegLat;
  const double dx = bx - ax, dy = by - ay;
  const double len2 = dx * dx + dy * dy;
  double t = len2 > 0 ? -(ax * dx + ay * dy) / len2 : 0;
  t = std::clamp(t, 0.0, 1.0);
  const double x = ax + t * dx, y = ay + t * dy;
  return std::sqrt(x * x + y * y);
}

// Lower bound, in meters, of the distance from p to anything inside box.
double BoxDistance(const PointLL& p, const AABB2<PointLL>& box) {
  const PointLL closest(std::clamp(p.lng(), box.minx(), box.maxx()),
                        std::clamp(p.lat(), box.miny(), box.maxy()));
  return Meters(p, closest) * kLowerBoundFactor;
}

std::string Trim(const std::string& s) {
  const auto first = s.find_first_not_of(" \t\r");
  if (first == std::string::npos) {
    return {};
  }
  return s.substr(first, s.find_last_not_of(" \t\r") - first + 1);
}

} // namespace

namespace valhalla {
namespace baldr {

// ---------------------------------------------------------------------------------------------
// RegionPolygon

int32_t RegionPolygon::Row(double lat) {
  return static_cast<int32_t>(std::floor(lat / kCellDeg));
}

int32_t RegionPolygon::Col(double lng) {
  return static_cast<int32_t>(std::floor(lng / kCellDeg));
}

RegionPolygon::RegionPolygon(const std::string& poly_file) {
  std::ifstream in(poly_file);
  if (!in) {
    throw std::runtime_error("Cannot open polygon file " + poly_file);
  }
  size_t line_number = 1;
  const auto fail = [&](const std::string& what) {
    return std::runtime_error("Malformed polygon file " + poly_file + " line " +
                              std::to_string(line_number) + ": " + what);
  };
  std::string line;
  if (!std::getline(in, line)) {
    throw fail("empty file");
  }
  // the first line is the polygon name; then rings (a name line, coordinate lines, END), then END
  std::vector<PointLL> ring;
  bool in_ring = false;
  auto close_ring = [&]() {
    if (ring.size() < 3) {
      throw fail("ring with fewer than three points");
    }
    for (size_t i = 0; i + 1 < ring.size(); ++i) {
      segments_.push_back({ring[i], ring[i + 1]});
    }
    if (!(ring.front() == ring.back())) {
      segments_.push_back({ring.back(), ring.front()});
    }
    ring.clear();
  };
  while (std::getline(in, line)) {
    ++line_number;
    line = Trim(line);
    if (line.empty()) {
      continue;
    }
    if (line == "END") {
      if (in_ring) {
        close_ring();
        in_ring = false;
      }
      continue;
    }
    if (!in_ring) {
      in_ring = true; // ring name, prefixed with '!' for holes
      continue;
    }
    std::istringstream ss(line);
    double lng, lat;
    if (!(ss >> lng >> lat)) {
      throw fail("expected a longitude and a latitude");
    }
    if (!std::isfinite(lng) || !std::isfinite(lat) || std::abs(lng) > 180 || std::abs(lat) > 90) {
      throw fail("coordinate out of range");
    }
    ring.emplace_back(lng, lat);
  }
  if (in_ring) {
    throw fail("ring without END");
  }
  if (segments_.empty()) {
    throw fail("no rings");
  }

  bbox_ = AABB2<PointLL>(segments_.front().a, segments_.front().a);
  for (uint32_t i = 0; i < segments_.size(); ++i) {
    const auto& s = segments_[i];
    bbox_.Expand(s.a);
    bbox_.Expand(s.b);
    const int32_t r0 = Row(std::min(s.a.lat(), s.b.lat())), r1 = Row(std::max(s.a.lat(), s.b.lat()));
    const int32_t c0 = Col(std::min(s.a.lng(), s.b.lng())), c1 = Col(std::max(s.a.lng(), s.b.lng()));
    row_min_ = std::min(row_min_, r0);
    row_max_ = std::max(row_max_, r1);
    col_min_ = std::min(col_min_, c0);
    col_max_ = std::max(col_max_, c1);
    for (int32_t r = r0; r <= r1; ++r) {
      rows_[r].push_back(i);
      for (int32_t c = c0; c <= c1; ++c) {
        cells_[Key(r, c)].push_back(i);
      }
    }
  }
}

bool RegionPolygon::Contains(const PointLL& p) const {
  if (!bbox_.Contains(p)) {
    return false;
  }
  const auto row = rows_.find(Row(p.lat()));
  if (row == rows_.end()) {
    return false;
  }
  // even-odd ray cast towards increasing longitude
  bool inside = false;
  for (auto i : row->second) {
    const auto& a = segments_[i].a;
    const auto& b = segments_[i].b;
    if ((a.lat() > p.lat()) != (b.lat() > p.lat())) {
      const double x = a.lng() + (p.lat() - a.lat()) * (b.lng() - a.lng()) / (b.lat() - a.lat());
      if (x > p.lng()) {
        inside = !inside;
      }
    }
  }
  return inside;
}

double RegionPolygon::Distance(const PointLL& p) const {
  const int32_t row = Row(p.lat()), col = Col(p.lng());
  // narrowest cell width near p, for the bound on segments outside the searched rings
  const double min_cell_m =
      kCellDeg *
      std::min(kMetersPerDegLat,
               kMetersPerDegLng * std::cos((std::abs(p.lat()) + kCellDeg * 8) * kPiD / 180.0));
  double best = std::numeric_limits<double>::max();
  auto visit = [&](int32_t r, int32_t c) {
    const auto cell = cells_.find(Key(r, c));
    if (cell != cells_.end()) {
      for (auto i : cell->second) {
        best = std::min(best, SegmentDistance(p, segments_[i].a, segments_[i].b));
      }
    }
  };
  // cells outside the range of the index hold no segment, so a ring is visited only where it
  // overlaps the range: far from the polygon, most of a ring is outside it
  const auto visit_row = [&](int32_t r, int32_t c0, int32_t c1) {
    if (r >= row_min_ && r <= row_max_) {
      for (int32_t c = std::max(c0, col_min_); c <= std::min(c1, col_max_); ++c) {
        visit(r, c);
      }
    }
  };
  const auto visit_col = [&](int32_t c, int32_t r0, int32_t r1) {
    if (c >= col_min_ && c <= col_max_) {
      for (int32_t r = std::max(r0, row_min_); r <= std::min(r1, row_max_); ++r) {
        visit(r, c);
      }
    }
  };
  for (int32_t ring = 0; ring <= kMaxDistanceRings; ++ring) {
    // only the ring's perimeter: its interior was visited by the smaller rings
    if (ring == 0) {
      visit_row(row, col, col);
    } else {
      visit_row(row - ring, col - ring, col + ring);
      visit_row(row + ring, col - ring, col + ring);
      visit_col(col - ring, row - ring + 1, row + ring - 1);
      visit_col(col + ring, row - ring + 1, row + ring - 1);
    }
    // everything outside this ring is at least ring cells away
    if (best <= ring * min_cell_m) {
      return best;
    }
  }
  for (const auto& s : segments_) {
    best = std::min(best, SegmentDistance(p, s.a, s.b));
  }
  return best;
}

// ---------------------------------------------------------------------------------------------
// PackageSet

namespace {

// FNV-1a: a hash that is the same on every platform, for join keys and tar fingerprints.
class Hasher {
public:
  void Bytes(const void* data, size_t size) {
    const auto* bytes = static_cast<const unsigned char*>(data);
    for (size_t i = 0; i < size; ++i) {
      hash_ = (hash_ ^ bytes[i]) * 1099511628211ull;
    }
  }
  void U64(uint64_t value) {
    unsigned char bytes[8];
    for (int i = 0; i < 8; ++i) {
      bytes[i] = static_cast<unsigned char>(value >> (8 * i));
    }
    Bytes(bytes, sizeof(bytes));
  }
  void String(const std::string& value) {
    U64(value.size());
    Bytes(value.data(), value.size());
  }
  uint64_t value() const {
    return hash_;
  }

private:
  uint64_t hash_ = 14695981039346656037ull;
};

// Little-endian fields of the binary formats.
void PutU32(std::string& out, uint32_t value) {
  for (int i = 0; i < 4; ++i) {
    out.push_back(static_cast<char>(value >> (8 * i)));
  }
}
void PutU64(std::string& out, uint64_t value) {
  for (int i = 0; i < 8; ++i) {
    out.push_back(static_cast<char>(value >> (8 * i)));
  }
}
uint64_t GetLE(const char* data, int bytes) {
  uint64_t value = 0;
  for (int i = bytes - 1; i >= 0; --i) {
    value = (value << 8) | static_cast<unsigned char>(data[i]);
  }
  return value;
}

uint64_t HashFile(const std::string& file) {
  std::ifstream in(file, std::ios::binary);
  if (!in) {
    throw std::runtime_error("Cannot open " + file);
  }
  Hasher hasher;
  char buffer[1 << 16];
  while (in.read(buffer, sizeof(buffer)) || in.gcount() > 0) {
    hasher.Bytes(buffer, static_cast<size_t>(in.gcount()));
  }
  return hasher.value();
}

// A package tar held open from the moment its read-only mapping is made, so the private mapping
// of its first rewrite maps the same file even when the path has been replaced since.
class OpenFile {
public:
  explicit OpenFile(const std::string& file) : file_(file) {
#ifdef _WIN32
    fd_ = _open(file.c_str(), _O_RDONLY | _O_BINARY);
#else
    fd_ = open(file.c_str(), O_RDONLY | O_CLOEXEC);
#endif
    if (fd_ == -1) {
      throw std::runtime_error(file + " (open): " + strerror(errno));
    }
    struct stat s;
    if (fstat(fd_, &s) != 0) {
      const int error = errno;
      Close();
      throw std::runtime_error(file + " (fstat): " + strerror(error));
    }
    size_ = static_cast<size_t>(s.st_size);
    device_ = static_cast<uint64_t>(s.st_dev);
    inode_ = static_cast<uint64_t>(s.st_ino);
  }
  ~OpenFile() {
    Close();
  }
  OpenFile(const OpenFile&) = delete;
  OpenFile& operator=(const OpenFile&) = delete;

  // Whether the path still names the open file, with the given size.
  bool Is(const std::string& path, size_t size) const {
    if (size != size_) {
      return false;
    }
#ifdef _WIN32
    (void)path;
    return true; // no inode numbers; files that are in use cannot be replaced there
#else
    struct stat s;
    return stat(path.c_str(), &s) == 0 && static_cast<uint64_t>(s.st_dev) == device_ &&
           static_cast<uint64_t>(s.st_ino) == inode_;
#endif
  }

  int fd() const {
    return fd_;
  }
  // Reads size bytes at offset.
  void Read(char* data, size_t size, size_t offset) const {
    while (size > 0) {
#ifdef _WIN32
      const auto n = _lseeki64(fd_, static_cast<__int64>(offset), SEEK_SET) < 0
                         ? -1
                         : _read(fd_, data, static_cast<unsigned>(std::min<size_t>(size, 1u << 30)));
#else
      const ssize_t n = pread(fd_, data, size, static_cast<off_t>(offset));
      if (n < 0 && errno == EINTR) {
        continue;
      }
#endif
      if (n <= 0) {
        throw std::runtime_error(file_ + " (read): " + (n < 0 ? strerror(errno) : "unexpected end"));
      }
      data += n;
      size -= static_cast<size_t>(n);
      offset += static_cast<size_t>(n);
    }
  }
  size_t size() const {
    return size_;
  }
  const std::string& file() const {
    return file_;
  }

private:
  void Close() {
#ifdef _WIN32
    _close(fd_);
#else
    close(fd_);
#endif
  }

  std::string file_;
  int fd_ = -1;
  size_t size_ = 0;
  uint64_t device_ = 0, inode_ = 0;
};

// A package tar mapped copy-on-write, for the runtime join: tiles are rewritten in place, and only
// the pages actually written become private memory. The file itself is never modified.
struct PrivateMap {
  explicit PrivateMap(const OpenFile& file) : size(file.size()) {
    if (size == 0) {
      throw std::runtime_error(file.file() + " is empty");
    }
    void* ptr = mmap(nullptr, size, PROT_READ | PROT_WRITE, MAP_PRIVATE, file.fd(), 0);
    if (ptr == MAP_FAILED) {
      throw std::runtime_error(file.file() + " (mmap): " + strerror(errno));
    }
    data = static_cast<char*>(ptr);
  }
  ~PrivateMap() {
    munmap(data, size);
  }
  PrivateMap(const PrivateMap&) = delete;
  PrivateMap& operator=(const PrivateMap&) = delete;

  char* data = nullptr;
  size_t size = 0;
};

// Tile memory inside a mapping, which it keeps alive.
template <typename Mapping> class MappedTileMemory final : public GraphMemory {
public:
  MappedTileMemory(std::shared_ptr<Mapping> mapping, const char* tile, size_t tile_size)
      : mapping_(std::move(mapping)) {
    data = const_cast<char*>(tile);
    size = tile_size;
  }

private:
  std::shared_ptr<Mapping> mapping_;
};

// Tile ids of package copies beyond the first, process-wide, because GraphTile::BoundingBox and
// TileHierarchy::GetGraphIdBoundingBox map them back to real tiles without a PackageSet. Each
// package set registers the ids it uses and releases them when it is destroyed; an id can only be
// held for one (real tile, copy index) at a time.
//
// Lookup is lock-free, as it runs on every bounding box query; changes take the mutex. Entries
// are published with release stores after their chunk.
class CopyIdTable {
public:
  // Never destroyed: package sets held by static objects may unregister during exit.
  static CopyIdTable& Get() {
    static CopyIdTable* table = new CopyIdTable();
    return *table;
  }

  // (real tile id, copy index) of a tile id; (tileid, 0) for real and unregistered ids.
  std::pair<uint32_t, uint32_t> Lookup(uint32_t level, uint32_t tileid) const {
    const uint64_t entry = Entry(level, tileid);
    if (entry == 0) {
      return {tileid, 0};
    }
    return {static_cast<uint32_t>(entry), static_cast<uint32_t>(entry >> 32)};
  }

  // Encoded (real tile id, copy index) an id is registered for, or 0.
  uint64_t Entry(uint32_t level, uint32_t tileid) const {
    if (level >= levels_.size() || tileid < levels_[level].first || tileid > kMaxGraphTileId) {
      return 0;
    }
    const auto& l = levels_[level];
    const uint32_t index = tileid - l.first;
    const auto* entries = l.chunks[index >> kChunkBits].load(std::memory_order_acquire);
    return entries ? entries[index & (kChunkSize - 1)].load(std::memory_order_acquire) : 0;
  }

  // The caller holds mutex(). Takes one reference to the id for the entry; false if the id is
  // registered for another entry.
  bool Acquire(uint32_t level, uint32_t tileid, uint64_t entry) {
    auto& l = levels_.at(level);
    const uint32_t index = tileid - l.first;
    auto& chunk = l.chunks[index >> kChunkBits];
    auto* entries = chunk.load(std::memory_order_acquire);
    if (!entries) {
      entries = new std::atomic<uint64_t>[kChunkSize]();
      chunk.store(entries, std::memory_order_release);
    }
    auto& slot = entries[index & (kChunkSize - 1)];
    const uint64_t current = slot.load(std::memory_order_acquire);
    if (current != 0 && current != entry) {
      return false;
    }
    slot.store(entry, std::memory_order_release);
    ++l.references[tileid];
    return true;
  }

  // The caller holds mutex(). Drops one reference to the id, freeing it with the last one.
  void Release(uint32_t level, uint32_t tileid) {
    auto& l = levels_.at(level);
    const auto references = l.references.find(tileid);
    if (references == l.references.end() || --references->second > 0) {
      return;
    }
    l.references.erase(references);
    const uint32_t index = tileid - l.first;
    l.chunks[index >> kChunkBits]
        .load(std::memory_order_acquire)[index & (kChunkSize - 1)]
        .store(0, std::memory_order_release);
  }

  std::mutex& mutex() {
    return mutex_;
  }

private:
  static constexpr uint32_t kChunkBits = 12;
  static constexpr uint32_t kChunkSize = 1u << kChunkBits;

  struct Level {
    uint32_t first = 0; // first id above the real grid
    // id - first -> (copy << 32 | real tile id), in lazily allocated chunks
    std::unique_ptr<std::atomic<std::atomic<uint64_t>*>[]> chunks;
    std::unordered_map<uint32_t, uint32_t> references; // id -> package sets holding it
  };

  CopyIdTable() : levels_(TileHierarchy::levels().size()) {
    for (size_t level = 0; level < levels_.size(); ++level) {
      auto& l = levels_[level];
      l.first = TileHierarchy::levels()[level].tiles.TileCount();
      const uint32_t chunk_count = ((kMaxGraphTileId + 1 - l.first) >> kChunkBits) + 1;
      l.chunks.reset(new std::atomic<std::atomic<uint64_t>*>[chunk_count]());
    }
  }

  std::mutex mutex_;
  std::vector<Level> levels_;
};

// Copy ids of a level. Every copy of every tile up to FixedCopies has a fixed id, the same in
// every package set: copy c of tile t gets TileCount * c + t. Further copies get dynamic ids above
// those, per set. The largest tile id is never handed out: references into tiles a package does
// not ship point there.
uint32_t FixedCopies(uint32_t level) {
  const uint32_t count = TileHierarchy::levels()[level].tiles.TileCount();
  return (kMaxGraphTileId - count) / count;
}
uint32_t FirstDynamicId(uint32_t level) {
  return TileHierarchy::levels()[level].tiles.TileCount() * (FixedCopies(level) + 1);
}

// A GraphId in no tile of any package set, for references into tiles a package does not ship.
valhalla::baldr::GraphId Nowhere(const valhalla::baldr::GraphId& id) {
  return {valhalla::baldr::kMaxGraphTileId, id.level(), id.id()};
}

// Whether id is the base id of a tile of the (non transit) tile grid.
bool IsGridTile(const valhalla::baldr::GraphId& id) {
  return id.is_valid() && id.id() == 0 &&
         id.level() < valhalla::baldr::TileHierarchy::levels().size() &&
         id.tileid() < valhalla::baldr::TileHierarchy::levels()[id.level()].tiles.TileCount();
}

// Calls fn for every complex restriction in a tile's raw restriction array.
template <typename Fn> void ForEachRestriction(char* data, size_t size, Fn fn) {
  for (size_t offset = 0; offset < size;) {
    if (size - offset < sizeof(valhalla::baldr::ComplexRestriction)) {
      throw std::runtime_error("Truncated complex restriction");
    }
    auto* restriction = reinterpret_cast<valhalla::baldr::ComplexRestriction*>(data + offset);
    const size_t length = restriction->SizeOf();
    if (length > size - offset) {
      throw std::runtime_error("Truncated complex restriction");
    }
    fn(*restriction);
    offset += length;
  }
}

// Index of the node of a tile that an edge of the tile starts at: the last node whose first edge
// index is not after the edge. -1 if there is none.
int64_t StartNode(const valhalla::baldr::GraphTile& tile, uint32_t edge) {
  const auto* nodes = tile.node(0);
  const auto* node = std::upper_bound(nodes, nodes + tile.header()->nodecount(), edge,
                                      [](uint32_t e, const valhalla::baldr::NodeInfo& n) {
                                        return e < n.edge_index();
                                      });
  return node == nodes ? -1 : node - 1 - nodes;
}

// A tile inside a package tar.
struct TarTile {
  uint64_t id; // tile base
  size_t offset, size;
};

// Lists the tiles of a tar in id order and returns the tar's fingerprint: a hash of every tile's
// id, size, and header. Tile headers hold a hash of the tile's data (GraphTileHeader::
// tile_checksum), so the fingerprint changes with the contents of any tile without reading them;
// a tile whose header holds no hash is hashed whole.
uint64_t ScanTar(midgard::tar& archive, const std::string& file, std::vector<TarTile>& tiles) {
  using valhalla::baldr::GraphId;
  using valhalla::baldr::GraphTileHeader;
  const char* base = archive.mm.get();
  archive.for_each([&](const std::string& name, const char* data, size_t size) {
    if (std::filesystem::path(name).extension() != ".gph") {
      return true; // not a tile
    }
    GraphId id;
    try {
      id = GraphId::FromTilePath(name);
    } catch (const std::exception&) { return true; }
    if (!IsGridTile(id)) {
      throw std::runtime_error(file + " has tile " + name + ": only levels 0 to " +
                               std::to_string(valhalla::baldr::TileHierarchy::levels().size() - 1) +
                               " of the tile grid are supported");
    }
    const auto offset = static_cast<size_t>(data - base);
    if (size < sizeof(GraphTileHeader) || size > archive.mm.size() - offset) {
      throw std::runtime_error(file + " has a truncated tile " + name);
    }
    tiles.push_back({id.value, offset, size});
    return true;
  });
  if (tiles.empty()) {
    throw std::runtime_error(file + " has no tiles");
  }
  std::sort(tiles.begin(), tiles.end(),
            [](const TarTile& a, const TarTile& b) { return a.id < b.id; });
  Hasher hasher;
  hasher.U64(tiles.size());
  for (size_t i = 0; i < tiles.size(); ++i) {
    if (i > 0 && tiles[i].id == tiles[i - 1].id) {
      throw std::runtime_error(file + " has tile " + std::to_string(GraphId(tiles[i].id)) + " twice");
    }
    hasher.U64(tiles[i].id);
    hasher.U64(tiles[i].size);
    const char* tile = base + tiles[i].offset;
    // tiles written by other tile builders may lack the data hash: hash their data instead
    GraphTileHeader header;
    std::memcpy(&header, tile, sizeof(header));
    hasher.Bytes(tile, header.tile_checksum() != 0 ? sizeof(GraphTileHeader) : tiles[i].size);
  }
  return hasher.value();
}

// Bake-time references of one tile (see PackageSet::WriteTileRefs).
struct TileRefs {
  AABB2<PointLL> reach; // nodes outside the tile its GraphIds lead to (if has_reach)
  bool has_reach = false;
  std::vector<uint64_t> ids; // other tiles its GraphIds point into
};

// Reads a tile reference index. Little endian: uint32 magic, uint32 format version, uint64
// fingerprint of the indexed tar, uint32 tile count, then per tile a uint32 tile base id, a float
// reach box (min lng, min lat, max lng, max lat; min > max when empty), a uint32 reference count
// and that many uint32 tile base ids. An index that cannot be read or does not match the package
// is ignored with a warning: every tile of the package is then joined in full.
std::unordered_map<uint64_t, TileRefs> ReadTileRefs(const std::string& file,
                                                    const std::string& package,
                                                    uint64_t fingerprint,
                                                    const std::unordered_set<uint64_t>& tiles) {
  using valhalla::baldr::GraphId;
  std::unordered_map<uint64_t, TileRefs> index;
  try {
    std::ifstream in(file, std::ios::binary);
    if (!in) {
      throw std::runtime_error("cannot open it");
    }
    const auto read = [&](void* target, size_t size) {
      if (!in.read(static_cast<char*>(target), static_cast<std::streamsize>(size))) {
        throw std::runtime_error("truncated");
      }
    };
    const auto read_u32 = [&]() {
      char bytes[4];
      read(bytes, sizeof(bytes));
      return static_cast<uint32_t>(GetLE(bytes, 4));
    };
    const auto read_tile = [&]() {
      const GraphId id(read_u32());
      if (!IsGridTile(id)) {
        throw std::runtime_error("invalid tile id " + std::to_string(id.value));
      }
      return id.value;
    };
    if (read_u32() != kTileRefsMagic) {
      throw std::runtime_error("not a tile reference index");
    }
    if (const auto version = read_u32(); version != kTileRefsVersion) {
      throw std::runtime_error("format version " + std::to_string(version) + ", expected " +
                               std::to_string(kTileRefsVersion));
    }
    char stored[8];
    read(stored, sizeof(stored));
    if (GetLE(stored, 8) != fingerprint) {
      throw std::runtime_error("built for another tile tar");
    }
    const uint32_t count = read_u32();
    index.reserve(count);
    for (uint32_t t = 0; t < count; ++t) {
      const auto [refs, inserted] = index.emplace(read_tile(), TileRefs{});
      if (!inserted) {
        throw std::runtime_error("duplicate tile " + std::to_string(GraphId(refs->first)));
      }
      if (!tiles.count(refs->first)) {
        throw std::runtime_error("lists tile " + std::to_string(GraphId(refs->first)) +
                                 " that the package does not have");
      }
      float box[4];
      for (auto& value : box) {
        const uint32_t bits = read_u32();
        std::memcpy(&value, &bits, sizeof(value));
      }
      refs->second.has_reach = box[0] <= box[2];
      if (refs->second.has_reach) {
        refs->second.reach = AABB2<PointLL>(box[0], box[1], box[2], box[3]);
      }
      const uint32_t ids = read_u32();
      refs->second.ids.reserve(ids);
      for (uint32_t i = 0; i < ids; ++i) {
        refs->second.ids.push_back(read_tile());
      }
    }
    if (in.peek() != std::char_traits<char>::eof()) {
      throw std::runtime_error("trailing data");
    }
  } catch (const std::exception& e) {
    LOG_WARN("Ignoring tile reference index {} of package {} ({}): every tile of the package is "
             "joined in full",
             file, package, e.what());
    index.clear();
  }
  return index;
}

// Assigns only when the value changes: the tile may live in a private (copy-on-write) mapping,
// where every write, even of an identical value, turns a shared file page into private memory.
template <typename T> void Store(T& target, const T& value) {
  if (std::memcmp(&target, &value, sizeof(T)) != 0) {
    target = value;
  }
}

// Writes a file, durably once Sync is called.
class DurableFile {
public:
  explicit DurableFile(const std::string& file) : file_(file), out_(std::fopen(file.c_str(), "wb")) {
    if (!out_) {
      throw std::runtime_error("Cannot write " + file + ": " + strerror(errno));
    }
  }
  ~DurableFile() {
    if (out_) {
      std::fclose(out_);
    }
  }
  DurableFile(const DurableFile&) = delete;
  DurableFile& operator=(const DurableFile&) = delete;

  void Write(const std::string& bytes) {
    if (!bytes.empty() && std::fwrite(bytes.data(), 1, bytes.size(), out_) != bytes.size()) {
      throw std::runtime_error("Cannot write " + file_ + ": " + strerror(errno));
    }
  }

  void Sync() {
    bool ok = std::fflush(out_) == 0;
#ifdef _WIN32
    ok = ok && _commit(_fileno(out_)) == 0;
#else
    ok = ok && fsync(fileno(out_)) == 0;
#endif
    ok = std::fclose(out_) == 0 && ok;
    out_ = nullptr;
    if (!ok) {
      throw std::runtime_error("Cannot write " + file_ + ": " + strerror(errno));
    }
  }

private:
  std::string file_;
  std::FILE* out_;
};

// Makes renames in a directory durable; best effort.
void SyncDirectory(const std::string& dir) {
#ifndef _WIN32
  const int fd = open(dir.c_str(), O_RDONLY | O_CLOEXEC);
  if (fd != -1) {
    fsync(fd);
    close(fd);
  }
#else
  (void)dir;
#endif
}

} // namespace

GraphId PackageSet::Real(const GraphId& id) {
  if (!id.is_valid()) {
    return id;
  }
  return GraphId(CopyIdTable::Get().Lookup(id.level(), id.tileid()).first, id.level(), id.id());
}

uint32_t PackageSet::CopyIndex(const GraphId& id) {
  return id.is_valid() ? CopyIdTable::Get().Lookup(id.level(), id.tileid()).second : 0;
}

bool PackageSet::RegisterCopyIds(std::vector<CopyId>& ids) {
  auto& table = CopyIdTable::Get();
  std::lock_guard<std::mutex> lock(table.mutex());
  const auto key = [](uint32_t level, uint32_t id) { return (uint64_t(level) << 32) | id; };
  const auto entry = [](const CopyId& c) { return (uint64_t(c.copy) << 32) | c.tileid; };
  std::unordered_set<uint64_t> taken;
  for (const auto& c : ids) {
    taken.insert(key(c.level, c.id));
  }
  bool kept = true;
  for (auto& c : ids) {
    const uint64_t current = table.Entry(c.level, c.id);
    if (current == 0 || current == entry(c)) {
      continue;
    }
    if (c.id < FirstDynamicId(c.level)) {
      throw std::logic_error("Fixed copy id " + std::to_string(c.id) + " taken on level " +
                             std::to_string(c.level));
    }
    // another live set holds this dynamic id for another tile: take the lowest free one
    taken.erase(key(c.level, c.id));
    uint32_t id = FirstDynamicId(c.level);
    while (id < kMaxGraphTileId && (table.Entry(c.level, id) != 0 || taken.count(key(c.level, id)))) {
      ++id;
    }
    if (id == kMaxGraphTileId) {
      throw std::runtime_error("No free tile ids left on level " + std::to_string(c.level));
    }
    c.id = id;
    taken.insert(key(c.level, c.id));
    kept = false;
  }
  for (const auto& c : ids) {
    if (!table.Acquire(c.level, c.id, entry(c))) {
      throw std::logic_error("Copy id " + std::to_string(c.id) + " taken on level " +
                             std::to_string(c.level));
    }
  }
  return kept;
}

void PackageSet::UnregisterCopyIds(const std::vector<CopyId>& ids) {
  auto& table = CopyIdTable::Get();
  std::lock_guard<std::mutex> lock(table.mutex());
  for (const auto& c : ids) {
    table.Release(c.level, c.id);
  }
}

struct PackageSet::TileEntry {
  enum class State { kPending, kClean, kRewritten, kBroken };

  uint64_t real = 0;             // real tile base
  size_t offset = 0, size = 0;   // position in the package tar
  uint32_t copy = 0, tileid = 0; // copy index and tile id of this package's copy
  // position of the tile's patch in the overlay file; size 0 when the join left the tile as it is
  size_t overlay_offset = 0, overlay_size = 0;

  // the package's own tile, over the read-only mapping of its tar
  mutable std::once_flag raw_once;
  mutable graph_tile_ptr raw;
  // grid of node positions for twin lookups: node indexes, sorted by (cell key, index)
  mutable std::once_flag grid_once;
  mutable double grid_cell = 0;
  mutable std::vector<uint32_t> grid;
  // the tile as patched from the overlay
  mutable std::once_flag overlay_once;
  mutable graph_tile_ptr overlay_tile;
  mutable std::atomic<bool> counted{false};
  // the runtime join, rewriting the tile in the package's private mapping once
  mutable std::mutex mutex;
  mutable State state = State::kPending;
  mutable graph_tile_ptr joined;
};

struct PackageSet::Package {
  std::string name;
  int64_t build_time = 0;
  std::unique_ptr<RegionPolygon> polygon;
  uint64_t polygon_hash = 0;
  uint64_t fingerprint = 0; // of the tile tar
  std::string tile_extract;
  std::shared_ptr<midgard::tar> tar; // read-only mapping
  std::unique_ptr<OpenFile> file;    // the mapped tar, for the private mapping of rewritten tiles
  size_t tile_count = 0;
  std::unique_ptr<TileEntry[]> tiles;               // in real tile id order
  std::unordered_map<uint64_t, uint32_t> index;     // real tile base -> index in tiles
  std::unordered_map<uint64_t, TileRefs> refs;      // bake-time tile reference index; may be empty
  std::shared_ptr<const std::vector<char>> overlay; // when the set serves overlays
  mutable std::once_flag private_once;              // the private mapping, made on first rewrite
  mutable std::shared_ptr<PrivateMap> private_map;
  mutable std::mutex warned_mutex; // tiles the package references but does not ship, warned of
  mutable std::unordered_set<uint64_t> warned;

  const TileEntry* Find(uint64_t real) const {
    const auto found = index.find(real);
    return found == index.end() ? nullptr : &tiles[found->second];
  }
};

PackageSet::PackageSet() = default;

PackageSet::~PackageSet() {
  if (!registered_.empty()) {
    UnregisterCopyIds(registered_);
  }
}

std::unique_ptr<PackageSet> PackageSet::FromConfig(const boost::property_tree::ptree& pt) {
  const auto packages = pt.get_child_optional("packages");
  if (!packages || packages->empty()) {
    return nullptr;
  }
#ifndef ENABLE_THREAD_SAFE_TILE_REF_COUNT
  throw std::runtime_error("mjolnir.packages needs a build with ENABLE_THREAD_SAFE_TILE_REF_COUNT="
                           "ON: the tiles of a package set are shared by every thread");
#else
  auto set = Build(pt);
  set->Register();
  set->IndexCopies();
  const auto overlays = pt.get<std::string>("package_joined", "");
  if (!overlays.empty()) {
    set->OpenOverlays(overlays);
  }
  return set;
#endif
}

std::unique_ptr<PackageSet> PackageSet::Build(const boost::property_tree::ptree& pt) {
  std::unique_ptr<PackageSet> set(new PackageSet());
  set->join_tolerance_m_ = pt.get<double>("package_join_tolerance", 0.0);
  if (!(set->join_tolerance_m_ >= 0) || !std::isfinite(set->join_tolerance_m_)) {
    throw std::runtime_error("package_join_tolerance must be a non-negative number");
  }

  // in name order, so nothing depends on the order of the configuration
  std::vector<std::pair<std::string, const boost::property_tree::ptree*>> configs;
  for (const auto& kv : pt.get_child("packages")) {
    const auto name = kv.second.get<std::string>("name");
    // the name also names the package's overlay file
    if (name.empty() || name.find_first_of("/\\") != std::string::npos || name == "." ||
        name == "..") {
      throw std::runtime_error("Invalid package name '" + name + "'");
    }
    configs.emplace_back(name, &kv.second);
  }
  std::sort(configs.begin(), configs.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
  for (size_t i = 1; i < configs.size(); ++i) {
    if (configs[i].first == configs[i - 1].first) {
      throw std::runtime_error("Duplicate package name " + configs[i].first);
    }
  }

  for (const auto& [name, cfg] : configs) {
    auto pkg = std::make_unique<Package>();
    pkg->name = name;
    pkg->build_time = cfg->get<int64_t>("build_time", 0);
    const auto polygon = cfg->get<std::string>("polygon");
    pkg->polygon = std::make_unique<RegionPolygon>(polygon);
    pkg->polygon_hash = HashFile(polygon);
    pkg->tile_extract = cfg->get<std::string>("tile_extract");
    try {
      pkg->file = std::make_unique<OpenFile>(pkg->tile_extract);
      pkg->tar = std::make_shared<midgard::tar>(pkg->tile_extract);
    } catch (const std::exception& e) {
      throw std::runtime_error("Package " + name + ": " + e.what());
    }
    if (!pkg->file->Is(pkg->tile_extract, pkg->tar->mm.size())) {
      throw std::runtime_error("Package " + name + ": " + pkg->tile_extract +
                               " was replaced while it was loaded");
    }
    std::vector<TarTile> tiles;
    pkg->fingerprint = ScanTar(*pkg->tar, pkg->tile_extract, tiles);
    pkg->tile_count = tiles.size();
    pkg->tiles.reset(new TileEntry[tiles.size()]);
    std::unordered_set<uint64_t> ids;
    for (size_t i = 0; i < tiles.size(); ++i) {
      auto& entry = pkg->tiles[i];
      entry.real = tiles[i].id;
      entry.offset = tiles[i].offset;
      entry.size = tiles[i].size;
      pkg->index.emplace(entry.real, static_cast<uint32_t>(i));
      ids.insert(entry.real);
    }
    if (const auto refs = cfg->get_optional<std::string>("tile_refs")) {
      pkg->refs = ReadTileRefs(*refs, name, pkg->fingerprint, ids);
    }
    LOG_INFO("Package {} loaded: {} tiles, build_time {}", name, pkg->tile_count, pkg->build_time);
    set->packages_.push_back(std::move(pkg));
  }

  for (uint32_t p = 0; p < set->packages_.size(); ++p) {
    const auto& package = *set->packages_[p];
    for (size_t i = 0; i < package.tile_count; ++i) {
      set->copies_[package.tiles[i].real].push_back(p);
    }
  }

  // Slot 0 of a shared tile goes to the package owning the tile's center: routing there mostly
  // reads that copy, and references into slot 0 need no rewrite. Further copies follow in name
  // order. Tiles are visited in id order, so copy ids depend on nothing but the set.
  std::vector<uint64_t> tiles;
  tiles.reserve(set->copies_.size());
  for (const auto& kv : set->copies_) {
    tiles.push_back(kv.first);
  }
  std::sort(tiles.begin(), tiles.end());
  std::vector<uint32_t> next_dynamic(TileHierarchy::levels().size());
  for (uint32_t level = 0; level < next_dynamic.size(); ++level) {
    next_dynamic[level] = FirstDynamicId(level);
  }
  for (const auto tile_base : tiles) {
    auto& holders = set->copies_[tile_base];
    const GraphId id(tile_base);
    if (holders.size() > 1) {
      const auto center = TileHierarchy::GetGraphIdBoundingBox(id).Center();
      const auto it = std::find(holders.begin(), holders.end(), set->FindOwner(center));
      if (it != holders.end()) {
        std::rotate(holders.begin(), it, it + 1);
      }
    }
    const uint32_t count = TileHierarchy::levels()[id.level()].tiles.TileCount();
    for (uint32_t copy = 0; copy < holders.size(); ++copy) {
      uint32_t tileid = id.tileid();
      if (copy > 0 && copy <= FixedCopies(id.level())) {
        tileid = count * copy + id.tileid();
      } else if (copy > 0) {
        tileid = next_dynamic[id.level()]++;
        if (tileid >= kMaxGraphTileId) {
          throw std::runtime_error("No free tile ids left on level " + std::to_string(id.level()));
        }
      }
      auto& package = *set->packages_[holders[copy]];
      auto& entry = package.tiles[package.index.at(tile_base)];
      entry.copy = copy;
      entry.tileid = tileid;
    }
  }

  Hasher key;
  key.String("valhalla package join");
  key.U64(kJoinVersion);
  key.String(VALHALLA_VERSION);
  uint64_t tolerance_bits = 0;
  std::memcpy(&tolerance_bits, &set->join_tolerance_m_, sizeof(tolerance_bits));
  key.U64(tolerance_bits);
  key.U64(set->packages_.size());
  for (const auto& package : set->packages_) {
    key.String(package->name);
    key.U64(static_cast<uint64_t>(package->build_time));
    key.U64(package->polygon_hash);
    key.U64(package->fingerprint);
  }
  set->join_key_ = key.value();
  return set;
}

std::vector<PackageSet::CopyId> PackageSet::CopyIds() const {
  std::vector<CopyId> ids;
  for (const auto& package : packages_) {
    for (size_t i = 0; i < package->tile_count; ++i) {
      const auto& entry = package->tiles[i];
      if (entry.copy > 0) {
        const GraphId real(entry.real);
        ids.push_back({real.level(), real.tileid(), entry.copy, entry.tileid});
      }
    }
  }
  std::sort(ids.begin(), ids.end(), [](const CopyId& a, const CopyId& b) {
    return std::tie(a.level, a.tileid, a.copy) < std::tie(b.level, b.tileid, b.copy);
  });
  return ids;
}

void PackageSet::Register() {
  auto ids = CopyIds();
  if (!RegisterCopyIds(ids)) {
    LOG_WARN("Another package set in this process holds tile ids this one would use; it gets "
             "other ids and joins its tiles at runtime");
    std::unordered_map<uint64_t, uint32_t> replaced;
    for (const auto& c : ids) {
      replaced[GraphId(c.tileid, c.level, 0).value | (uint64_t(c.copy) << 46)] = c.id;
    }
    for (auto& package : packages_) {
      for (size_t i = 0; i < package->tile_count; ++i) {
        auto& entry = package->tiles[i];
        if (entry.copy > 0) {
          entry.tileid = replaced.at(entry.real | (uint64_t(entry.copy) << 46));
        }
      }
    }
  }
  registered_ = std::move(ids);
}

void PackageSet::IndexCopies() {
  by_copy_.clear();
  for (uint32_t p = 0; p < packages_.size(); ++p) {
    const auto& package = *packages_[p];
    for (uint32_t i = 0; i < package.tile_count; ++i) {
      const auto& entry = package.tiles[i];
      by_copy_[GraphId(entry.tileid, GraphId(entry.real).level(), 0).value] = {p, i};
    }
  }
}

void PackageSet::OpenOverlays(const std::string& dir) {
  const auto ids = CopyIds();
  // (package, tile index, offset, size) of every patch, applied once all overlays match
  std::vector<std::tuple<uint32_t, uint32_t, size_t, size_t>> tiles;
  std::vector<std::shared_ptr<const std::vector<char>>> overlays;
  const auto stale = [&](const std::string& file, const std::string& why) {
    LOG_WARN("Overlay {} {}: joining packages at runtime", file, why);
  };
  for (uint32_t p = 0; p < packages_.size(); ++p) {
    const auto& package = *packages_[p];
    const auto file = (std::filesystem::path(dir) / (package.name + kOverlaySuffix)).string();
    std::shared_ptr<std::vector<char>> overlay;
    try {
      if (!std::filesystem::exists(file)) {
        stale(file, "is missing");
        return;
      }
      std::ifstream in(file, std::ios::binary);
      overlay = std::make_shared<std::vector<char>>(std::filesystem::file_size(file));
      if (!in.read(overlay->data(), static_cast<std::streamsize>(overlay->size()))) {
        throw std::runtime_error(strerror(errno));
      }
    } catch (const std::exception& e) {
      stale(file, std::string("cannot be read (") + e.what() + ")");
      return;
    }
    const char* data = overlay->data();
    const size_t size = overlay->size();
    std::string why;
    uint64_t tile_count = 0;
    if (size < kOverlayHeaderSize || GetLE(data, 4) != kOverlayMagic ||
        GetLE(data + 4, 4) != kOverlayVersion) {
      why = "is not an overlay of this format";
    } else if (GetLE(data + 8, 8) != join_key_) {
      why = "was written for other packages";
    } else if (GetLE(data + kOverlayIdCountOffset, 4) != ids.size()) {
      why = "was written with other copy ids";
    } else {
      tile_count = GetLE(data + kOverlayTileCountOffset, 4);
      const size_t directory = kOverlayHeaderSize + kOverlayIdSize * ids.size();
      if (size < directory || tile_count > (size - directory) / kOverlayDirectorySize) {
        why = "is truncated";
      } else {
        Hasher hasher;
        hasher.Bytes(data + kOverlayHeaderSize, size - kOverlayHeaderSize);
        if (hasher.value() != GetLE(data + kOverlayHashOffset, 8)) {
          why = "is damaged";
        }
      }
      for (size_t i = 0; i < ids.size() && why.empty(); ++i) {
        const char* id = data + kOverlayHeaderSize + kOverlayIdSize * i;
        const CopyId stored{static_cast<uint32_t>(GetLE(id, 4)),
                            static_cast<uint32_t>(GetLE(id + 4, 4)),
                            static_cast<uint32_t>(GetLE(id + 8, 4)),
                            static_cast<uint32_t>(GetLE(id + 12, 4))};
        if (!(stored == ids[i])) {
          why = "was written with other copy ids";
        }
      }
    }
    if (why.empty()) {
      const size_t directory = kOverlayHeaderSize + kOverlayIdSize * ids.size();
      const size_t patches = directory + kOverlayDirectorySize * tile_count;
      uint64_t previous = 0, expected_offset = 0;
      for (uint64_t t = 0; t < tile_count && why.empty(); ++t) {
        const char* at = data + directory + kOverlayDirectorySize * t;
        const GraphId real(GetLE(at, 4));
        const size_t offset = GetLE(at + 4, 4), patch_size = GetLE(at + 8, 4);
        const auto found = package.index.find(real.value);
        if (t > 0 && real.value <= previous) {
          why = "has its tiles out of order";
        } else if (found == package.index.end()) {
          why = "has a tile the package does not have: " + std::to_string(real);
        } else if (offset != expected_offset || patch_size > size - patches - offset ||
                   patch_size < 8) {
          why = "is truncated";
        } else if (GetLE(data + patches + offset, 4) != package.tiles[found->second].size) {
          why = "has a patch for a tile of another size: " + std::to_string(real);
        } else {
          tiles.emplace_back(p, found->second, patches + offset, patch_size);
          previous = real.value;
          expected_offset = offset + patch_size;
        }
      }
      if (why.empty() && patches + expected_offset != size) {
        why = "has trailing bytes";
      }
    }
    if (!why.empty()) {
      stale(file, why);
      return;
    }
    overlays.push_back(std::move(overlay));
  }
  for (const auto& [p, i, offset, size] : tiles) {
    auto& entry = packages_[p]->tiles[i];
    entry.overlay_offset = offset;
    entry.overlay_size = size;
  }
  for (uint32_t p = 0; p < packages_.size(); ++p) {
    packages_[p]->overlay = std::move(overlays[p]);
  }
  overlays_ = true;
  LOG_INFO("Serving {} packages from the overlays in {}", packages_.size(), dir);
}

GraphId PackageSet::Virtual(uint32_t pkg, const GraphId& real) const {
  if (!real.is_valid()) {
    return real;
  }
  const auto* entry = packages_[pkg]->Find(real.tile_base().value);
  if (!entry) {
    WarnMissingTile(pkg, real.tile_base());
    return Nowhere(real);
  }
  return GraphId(entry->tileid, real.level(), real.id());
}

void PackageSet::WarnMissingTile(uint32_t pkg, const GraphId& tile) const {
  const auto& package = *packages_[pkg];
  std::lock_guard<std::mutex> lock(package.warned_mutex);
  if (package.warned.insert(tile.value).second) {
    LOG_WARN("Package {} references tile {} it does not ship: those references are disabled",
             package.name, std::to_string(tile));
  }
}

std::vector<GraphId> PackageSet::Copies(const GraphId& tile) const {
  std::vector<GraphId> result;
  uint64_t real = tile.tile_base().value;
  if (const auto copy = by_copy_.find(real); copy != by_copy_.end()) {
    real = packages_[copy->second.first]->tiles[copy->second.second].real;
  }
  const auto it = copies_.find(real);
  if (it != copies_.end()) {
    for (const uint32_t pkg : it->second) {
      result.emplace_back(packages_[pkg]->Find(real)->tileid, GraphId(real).level(), 0);
    }
  }
  return result;
}

std::vector<GraphId> PackageSet::AllTiles() const {
  std::vector<GraphId> result;
  result.reserve(by_copy_.size());
  for (const auto& kv : by_copy_) {
    result.emplace_back(kv.first);
  }
  std::sort(result.begin(), result.end());
  return result;
}

bool PackageSet::Exists(const GraphId& tile_base) const {
  return by_copy_.count(tile_base.tile_base().value) > 0;
}

int PackageSet::PackageOf(const GraphId& tile) const {
  const auto copy = by_copy_.find(tile.tile_base().value);
  return copy == by_copy_.end() ? -1 : static_cast<int>(copy->second.first);
}

size_t PackageSet::size() const {
  return packages_.size();
}

const std::string& PackageSet::name(size_t package) const {
  return packages_.at(package)->name;
}

PackageSet::Stats PackageSet::stats() const {
  Stats stats;
  stats.joined = counters_.joined;
  stats.lost = counters_.lost;
  stats.disabled = counters_.disabled;
  stats.clean_tiles = counters_.clean_tiles;
  stats.id_rewrites = counters_.id_rewrites;
  stats.full_rewrites = counters_.full_rewrites;
  stats.overlay_tiles = counters_.overlay_tiles;
  return stats;
}

int PackageSet::Owner(const PointLL& p) const {
  return FindOwner(p);
}

int PackageSet::FindOwner(const PointLL& p) const {
  // Depth is 1-Lipschitz, so when the depths at a level-2 tile's center differ by more than the
  // tile's diagonal, every point of the tile has the same owner.
  const auto& tiles = TileHierarchy::levels().back().tiles;
  const int32_t tile_id = tiles.TileId(p);
  if (tile_id < 0) {
    return PointOwner(p);
  }
  const int owner = CachedTileOwner(tile_id);
  if (owner != kUndecided) {
    return owner;
  }
  // undecided tile: try the ~1 km cell around p before computing the point itself
  const uint64_t key = CellKey(p, kOwnerCellDeg);
  int cell_owner = kUndecided;
  bool cached = false;
  {
    std::shared_lock<std::shared_mutex> lock(owner_mutex_);
    const auto cell = cell_owner_.find(key);
    if (cell != cell_owner_.end()) {
      cell_owner = cell->second;
      cached = true;
    }
  }
  if (!cached) {
    const auto row = std::floor(p.lat() / kOwnerCellDeg);
    const auto col = std::floor(p.lng() / kOwnerCellDeg);
    const AABB2<PointLL> box(col * kOwnerCellDeg, row * kOwnerCellDeg, (col + 1) * kOwnerCellDeg,
                             (row + 1) * kOwnerCellDeg);
    cell_owner = BoxOwner(box);
    std::unique_lock<std::shared_mutex> lock(owner_mutex_);
    cell_owner_.emplace(key, cell_owner);
  }
  return cell_owner == kUndecided ? PointOwner(p) : cell_owner;
}

int PackageSet::CachedTileOwner(int32_t tile_id) const {
  {
    std::shared_lock<std::shared_mutex> lock(owner_mutex_);
    const auto cached = tile_owner_.find(tile_id);
    if (cached != tile_owner_.end()) {
      return cached->second;
    }
  }
  // computed without the lock: the value depends only on the polygons
  const int owner = BoxOwner(TileHierarchy::levels().back().tiles.TileBounds(tile_id));
  std::unique_lock<std::shared_mutex> lock(owner_mutex_);
  tile_owner_.emplace(tile_id, owner);
  return owner;
}

bool PackageSet::NearOwnershipBoundary(const PointLL& p) const {
  const int32_t tile_id = TileHierarchy::levels().back().tiles.TileId(p);
  return tile_id < 0 || CachedTileOwner(tile_id) == kUndecided;
}

int PackageSet::BoxOwner(const AABB2<PointLL>& box) const {
  const PointLL center = box.Center();
  // half diagonal in meters, with slack for the local projection used by Distance()
  const double half_diag = kProjectionSlack * 0.5 *
                           Meters(PointLL(box.minx(), box.miny()), PointLL(box.maxx(), box.maxy()));
  // A single package whose boundary crosses the box owns all of it when every other package is
  // farther than it by more than the box diagonal: inside its polygon no one else is newer, and
  // outside it is the nearest package.
  int crossing = -1;
  for (uint32_t i = 0; i < packages_.size(); ++i) {
    const auto& polygon = *packages_[i]->polygon;
    if (polygon.bbox().Intersects(box) && polygon.Distance(center) <= half_diag) {
      if (crossing >= 0) {
        return kUndecided;
      }
      crossing = static_cast<int>(i);
    }
  }
  if (crossing >= 0) {
    const double d = packages_[crossing]->polygon->Distance(center);
    for (uint32_t i = 0; i < packages_.size(); ++i) {
      if (static_cast<int>(i) == crossing) {
        continue;
      }
      const auto& polygon = *packages_[i]->polygon;
      if (BoxDistance(center, polygon.bbox()) > d + 2 * half_diag + kDepthTieMeters) {
        continue;
      }
      if (polygon.Contains(center) ||
          polygon.Distance(center) <= d + 2 * half_diag + kDepthTieMeters) {
        return kUndecided;
      }
    }
    return crossing;
  }
  // No boundary crosses the box: every polygon contains all of it or none of it. The newest build
  // owns everything its polygon contains; depth only splits equal build times, and must do so for
  // the whole box.
  int best = -1, newest_count = 0;
  int64_t newest = 0;
  double best_depth = 0, second_depth = 0;
  for (uint32_t i = 0; i < packages_.size(); ++i) {
    const auto& polygon = *packages_[i]->polygon;
    if (!polygon.bbox().Intersects(box) || !polygon.Contains(center)) {
      continue;
    }
    const int64_t build_time = packages_[i]->build_time;
    if (best < 0 || build_time > newest) {
      best = static_cast<int>(i);
      newest = build_time;
      newest_count = 1;
      continue;
    }
    if (build_time < newest) {
      continue;
    }
    if (newest_count++ == 1) {
      best_depth = packages_[best]->polygon->Distance(center);
      second_depth = -1;
    }
    const double depth = polygon.Distance(center);
    if (depth > best_depth) {
      second_depth = best_depth;
      best_depth = depth;
      best = static_cast<int>(i);
    } else if (depth > second_depth) {
      second_depth = depth;
    }
  }
  if (best < 0) {
    // outside every polygon: the nearest package owns it, if that holds for the whole box
    double nearest_m = 0, second_m = 0;
    const int nearest = NearestPackage(center, nearest_m, second_m);
    return nearest >= 0 && nearest_m + 2 * half_diag + kDepthTieMeters < second_m ? nearest
                                                                                  : kUndecided;
  }
  if (newest_count == 1) {
    return best;
  }
  return best_depth - second_depth > 2 * half_diag + kDepthTieMeters ? best : kUndecided;
}

int PackageSet::NearestPackage(const PointLL& p, double& best, double& second) const {
  best = second = std::numeric_limits<double>::max();
  int nearest = -1;
  for (uint32_t i = 0; i < packages_.size(); ++i) {
    const auto& polygon = *packages_[i]->polygon;
    // the bounding box gives a lower bound; skip polygons that cannot come second
    if (BoxDistance(p, polygon.bbox()) > second) {
      continue;
    }
    const double d = polygon.Distance(p);
    if (d < best - kDepthTieMeters || (d <= best + kDepthTieMeters && nearest >= 0 &&
                                       packages_[i]->build_time > packages_[nearest]->build_time)) {
      second = std::min(second, best);
      best = d;
      nearest = static_cast<int>(i);
    } else {
      second = std::min(second, d);
    }
  }
  return nearest;
}

int PackageSet::PointOwner(const PointLL& p) const {
  // the newest build containing the point; among equal build times the polygon containing it most
  // deeply, then the lower index
  int best = -1;
  int64_t newest = 0;
  double best_depth = -1;
  for (uint32_t i = 0; i < packages_.size(); ++i) {
    if (!packages_[i]->polygon->Contains(p)) {
      continue;
    }
    const int64_t build_time = packages_[i]->build_time;
    if (best < 0 || build_time > newest) {
      best = static_cast<int>(i);
      newest = build_time;
      best_depth = -1;
      continue;
    }
    if (build_time < newest) {
      continue;
    }
    if (best_depth < 0) {
      best_depth = packages_[best]->polygon->Distance(p);
    }
    const double depth = packages_[i]->polygon->Distance(p);
    if (depth > best_depth + kDepthTieMeters) {
      best = static_cast<int>(i);
      best_depth = depth;
    }
  }
  if (best < 0) {
    // outside every package: the nearest package owns it, so every point has exactly one owner
    double nearest = 0, second = 0;
    return NearestPackage(p, nearest, second);
  }
  return best;
}

const GraphTile* PackageSet::RawTile(uint32_t pkg, const GraphId& real_base) const {
  const auto& package = *packages_[pkg];
  const auto* entry = package.Find(real_base.value);
  if (!entry) {
    return nullptr;
  }
  std::call_once(entry->raw_once, [&]() {
    const char* data = package.tar->mm.get() + entry->offset;
    entry->raw = GraphTile::Create(real_base,
                                   std::make_unique<MappedTileMemory<midgard::tar>>(package.tar, data,
                                                                                    entry->size));
  });
  return entry->raw.get();
}

bool PackageSet::EdgeOwned(uint32_t pkg, const GraphId& edge) const {
  const auto* tile = RawTile(pkg, edge.tile_base());
  if (!tile || edge.id() >= tile->header()->directededgecount()) {
    return false;
  }
  const int64_t node = StartNode(*tile, edge.id());
  return node >= 0 &&
         FindOwner(tile->node(static_cast<uint32_t>(node))->latlng(tile->header()->base_ll())) ==
             static_cast<int>(pkg);
}

bool PackageSet::NodeLL(uint32_t pkg, const GraphId& node, PointLL& ll) const {
  const auto* tile = RawTile(pkg, node.tile_base());
  if (!tile || node.id() >= tile->header()->nodecount()) {
    return false;
  }
  ll = tile->node(node.id())->latlng(tile->header()->base_ll());
  return true;
}

int64_t PackageSet::NearestNode(uint32_t pkg, const GraphId& tile_base, const PointLL& ll) const {
  const auto* tile = RawTile(pkg, tile_base);
  if (!tile) {
    return -1;
  }
  const auto& entry = *packages_[pkg]->Find(tile_base.value);
  const auto base_ll = tile->header()->base_ll();
  // keys are recomputed from the tile's nodes, so the grid costs 4 bytes per node
  const auto key_of = [&](uint32_t n) {
    return CellKey(tile->node(n)->latlng(base_ll), entry.grid_cell);
  };
  std::call_once(entry.grid_once, [&]() {
    // square cells at least as wide as the tolerance in both directions anywhere in (or just
    // outside) the tile, so a 3x3 search finds every node within the tolerance
    const auto box = TileHierarchy::GetGraphIdBoundingBox(tile_base);
    const double max_lat =
        std::min(90.0, std::max(std::abs(box.miny()), std::abs(box.maxy())) + kGridLatMarginDeg);
    const double meters_per_deg =
        std::min(kMetersPerDegLat, kMetersPerDegLng * std::cos(max_lat * kPiD / 180.0));
    entry.grid_cell = std::max(join_tolerance_m_, kMinGridCellMeters) /
                      std::max(meters_per_deg, kMinGridCellMeters);
    std::vector<std::pair<uint64_t, uint32_t>> keyed;
    keyed.reserve(tile->header()->nodecount());
    for (uint32_t n = 0; n < tile->header()->nodecount(); ++n) {
      keyed.emplace_back(key_of(n), n);
    }
    std::sort(keyed.begin(), keyed.end());
    entry.grid.reserve(keyed.size());
    for (const auto& k : keyed) {
      entry.grid.push_back(k.second);
    }
  });
  int64_t best = -1;
  double best_m = join_tolerance_m_ + 1e-6;
  for (int dr = -1; dr <= 1; ++dr) {
    for (int dc = -1; dc <= 1; ++dc) {
      const PointLL probe(ll.lng() + dc * entry.grid_cell, ll.lat() + dr * entry.grid_cell);
      const uint64_t key = CellKey(probe, entry.grid_cell);
      auto it = std::lower_bound(entry.grid.begin(), entry.grid.end(), key,
                                 [&](uint32_t n, uint64_t k) { return key_of(n) < k; });
      for (; it != entry.grid.end() && key_of(*it) == key; ++it) {
        const double m = Meters(tile->node(*it)->latlng(base_ll), ll);
        if (m < best_m || (m == best_m && *it < best)) {
          best_m = m;
          best = *it;
        }
      }
    }
  }
  return best;
}

int64_t
PackageSet::TwinNode(uint32_t from, const GraphId& node, const PointLL& ll, uint32_t to) const {
  const int64_t found = NearestNode(to, node.tile_base(), ll);
  if (found < 0) {
    return -1;
  }
  const auto* other = RawTile(to, node.tile_base());
  const PointLL found_ll =
      other->node(static_cast<uint32_t>(found))->latlng(other->header()->base_ll());
  return NearestNode(from, node.tile_base(), found_ll) == node.id() ? found : -1;
}

int PackageSet::MatchEdge(uint32_t to,
                          const GraphTile& tile,
                          uint32_t node,
                          uint64_t way_id,
                          uint32_t length,
                          uint32_t from,
                          const GraphId& far,
                          const PointLL& far_ll,
                          const char*& reason) const {
  // moved nodes change the length by up to twice the tolerance
  const double length_tolerance = std::max(2.0, 0.01 * length) + 2 * join_tolerance_m_;
  const NodeInfo* info = tile.node(node);
  reason = "no edge";
  for (uint32_t k = 0; k < info->edge_count(); ++k) {
    const DirectedEdge* edge = tile.directededge(info->edge_index() + k);
    if (edge->is_shortcut() || edge->IsTransitLine() || tile.edgeinfo(edge).wayid() != way_id) {
      continue;
    }
    PointLL end;
    if (!NodeLL(to, edge->endnode(), end) ||
        NearestNode(to, edge->endnode().tile_base(), far_ll) != edge->endnode().id() ||
        NearestNode(from, far.tile_base(), end) != far.id()) {
      reason = "span differs";
      continue;
    }
    if (std::abs(static_cast<double>(edge->length()) - length) > length_tolerance) {
      reason = "length differs";
      continue;
    }
    return static_cast<int>(k);
  }
  return -1;
}

bool PackageSet::FindTwin(uint32_t from_pkg,
                          const GraphTile& tile,
                          const DirectedEdge& edge,
                          const GraphId& start_tile,
                          uint32_t start_id,
                          const PointLL& start,
                          const PointLL& end,
                          uint32_t to_pkg,
                          Twin& twin) const {
  const GraphId end_real = edge.endnode();
  const auto* other = RawTile(to_pkg, end_real.tile_base());
  if (!other) {
    twin.reason = "no tile";
    return false;
  }
  // the other package's node nearest to our end node, and it must pick our end node back
  const int64_t found = TwinNode(from_pkg, end_real, end, to_pkg);
  if (found < 0) {
    twin.reason = "no node";
    return false;
  }
  // the edge back to our start node along the same way
  const int k =
      MatchEdge(to_pkg, *other, static_cast<uint32_t>(found), tile.edgeinfo(&edge).wayid(),
                edge.length(), from_pkg, GraphId(start_tile.tileid(), start_tile.level(), start_id),
                start, twin.reason);
  if (k < 0) {
    return false;
  }
  const DirectedEdge* back =
      other->directededge(other->node(static_cast<uint32_t>(found))->edge_index() + k);
  // the other package's copy of our edge is the opposing edge of `back`
  const auto* back_end_tile = RawTile(to_pkg, back->endnode().tile_base());
  const NodeInfo* back_end_node = back_end_tile->node(back->endnode().id());
  const DirectedEdge* same =
      back_end_tile->directededge(back_end_node->edge_index() + back->opp_index());

  twin.endnode =
      Virtual(to_pkg, GraphId(end_real.tileid(), end_real.level(), static_cast<uint32_t>(found)));
  twin.opp_index = static_cast<uint32_t>(k);
  twin.opp_local_idx = back->localedgeidx();
  twin.restrictions = same->restrictions();
  twin.reverseaccess = back->forwardaccess();
  return true;
}

GraphId PackageSet::RestrictionEdge(uint32_t pkg, const GraphId& real) const {
  if (!real.is_valid()) {
    return real;
  }
  const auto* tile = RawTile(pkg, real.tile_base());
  if (!tile) {
    WarnMissingTile(pkg, real.tile_base());
    return Nowhere(real);
  }
  const int64_t start =
      real.id() < tile->header()->directededgecount() ? StartNode(*tile, real.id()) : -1;
  if (start < 0) {
    LOG_WARN("Package {} restricts edge {} that it does not have", packages_[pkg]->name,
             std::to_string(real));
    return Nowhere(real);
  }
  const PointLL start_ll =
      tile->node(static_cast<uint32_t>(start))->latlng(tile->header()->base_ll());
  const int owner = FindOwner(start_ll);
  if (owner == static_cast<int>(pkg)) {
    return Virtual(pkg, real);
  }
  // routes use the owner's copy of the edge: the edge from the twin of its start node along the
  // same way to the twin of its end node
  const DirectedEdge* edge = tile->directededge(real.id());
  PointLL end_ll;
  const char* reason = "no node";
  if (NodeLL(pkg, edge->endnode(), end_ll)) {
    const GraphId start_node(real.tileid(), real.level(), static_cast<uint32_t>(start));
    const int64_t twin = TwinNode(pkg, start_node, start_ll, static_cast<uint32_t>(owner));
    if (twin >= 0) {
      const auto* other = RawTile(static_cast<uint32_t>(owner), real.tile_base());
      const int k = MatchEdge(static_cast<uint32_t>(owner), *other, static_cast<uint32_t>(twin),
                              tile->edgeinfo(edge).wayid(), edge->length(), pkg, edge->endnode(),
                              end_ll, reason);
      if (k >= 0) {
        const uint32_t index = other->node(static_cast<uint32_t>(twin))->edge_index() + k;
        return Virtual(static_cast<uint32_t>(owner), GraphId(real.tileid(), real.level(), index));
      }
    }
  }
  LOG_WARN("Package {}: no copy of restricted edge {} in package {} ({}), the restriction is not "
           "enforced",
           packages_[pkg]->name, std::to_string(real), packages_[owner]->name, reason);
  return Nowhere(real);
}

// The owners of the nodes of one tile, each found once: a node is asked for as the start of its
// edges and as the end of the edges that reach it, and the bins ask for the start nodes again.
class PackageSet::NodeOwners {
public:
  NodeOwners(const PackageSet& set, const GraphTile& tile)
      : set_(set), tile_(tile), owners_(tile.header()->nodecount(), kUnknown) {
  }

  int Of(uint32_t node) {
    int& owner = owners_[node];
    if (owner == kUnknown) {
      owner = set_.FindOwner(tile_.node(node)->latlng(tile_.header()->base_ll()));
    }
    return owner;
  }

private:
  static constexpr int kUnknown = std::numeric_limits<int>::min();
  const PackageSet& set_;
  const GraphTile& tile_;
  std::vector<int> owners_;
};

TileWords PackageSet::Layout(const GraphTile& tile) {
  static_assert(sizeof(NodeTransition) == 8 && sizeof(GraphId) == 8,
                "Tile patches count transitions and bin entries in 8-byte words");
  const auto* header = tile.header_;
  const char* base = reinterpret_cast<const char*>(header);
  // a region in words from the start of the tile; an unaligned region is left out of the patch's
  // classes, and its changes go into the general one
  const auto region = [&](const void* start, size_t count, size_t& first, size_t& words) {
    const auto offset = static_cast<size_t>(reinterpret_cast<const char*>(start) - base);
    if (offset % 8 == 0) {
      first = offset / 8;
      words = count;
    }
  };
  TileWords layout;
  region(tile.transitions_, header->transitioncount(), layout.transitions_first,
         layout.transitions_count);
  region(tile.directededges_, header->directededgecount(), layout.edges_first, layout.edges_count);
  region(tile.edge_bins_, header->bin_offset(kBinCount - 1).second, layout.bins_first,
         layout.bins_count);
  return layout;
}

void PackageSet::Rewrite(uint32_t pkg, GraphTile& tile, Work work) const {
  GraphTileHeader* header = tile.header_;
  const GraphId real_base = header->graphid();
  GraphTileHeader new_header = *header;
  new_header.set_graphid(Virtual(pkg, real_base));
  Store(*header, new_header);
  std::optional<NodeOwners> owners; // full rewrites only

  if (work == Work::kIdsOnly) {
    // every edge stays as it is, only end nodes in other tiles move to our slots
    for (uint32_t e = 0; e < header->directededgecount(); ++e) {
      DirectedEdge& target = tile.directededges_[e];
      DirectedEdge edge = target;
      edge.set_endnode(Virtual(pkg, edge.endnode()));
      Store(target, edge);
    }
  } else if (work == Work::kFull) {
    owners.emplace(*this, tile);
    JoinEdges(pkg, tile, real_base, *owners);
  }

  for (uint32_t t = 0; t < header->transitioncount(); ++t) {
    NodeTransition& transition = tile.transitions_[t];
    Store(transition, NodeTransition(Virtual(pkg, transition.endnode()), transition.up()));
  }

  for (size_t b = 0; b < kBinCount; ++b) {
    for (auto& id : tile.GetBin(b)) {
      if (!id.is_valid()) {
        continue;
      }
      if (work != Work::kFull) {
        // lighter work already knows every binned edge starts at a node of ours
        Store(id, Virtual(pkg, id));
      } else if (id.tile_base() == real_base) {
        // EdgeOwned, with the owners already found for this tile's nodes
        const int64_t node =
            id.id() < header->directededgecount() ? StartNode(tile, id.id()) : int64_t{-1};
        const bool owned = node >= 0 && owners->Of(static_cast<uint32_t>(node)) == static_cast<int>(pkg);
        Store(id, owned ? Virtual(pkg, id) : GraphId());
      } else if (!RawTile(pkg, id.tile_base())) {
        WarnMissingTile(pkg, id.tile_base());
        Store(id, GraphId());
      } else {
        Store(id, EdgeOwned(pkg, id) ? Virtual(pkg, id) : GraphId());
      }
    }
  }

  // lighter work already knows every restricted edge starts at a node of ours
  const auto edge_id = [&](const GraphId& id) {
    return work == Work::kFull ? RestrictionEdge(pkg, id) : Virtual(pkg, id);
  };
  const auto remap = [&](ComplexRestriction& restriction) {
    ComplexRestriction remapped = restriction;
    remapped.from_graphid_ = edge_id(restriction.from_graphid()).value;
    remapped.to_graphid_ = edge_id(restriction.to_graphid()).value;
    Store(restriction, remapped);
    auto* via = reinterpret_cast<GraphId*>(&restriction + 1);
    for (uint32_t v = 0; v < restriction.via_count(); ++v) {
      Store(via[v], edge_id(via[v]));
    }
  };
  ForEachRestriction(tile.complex_restriction_forward_, tile.complex_restriction_forward_size_,
                     remap);
  ForEachRestriction(tile.complex_restriction_reverse_, tile.complex_restriction_reverse_size_,
                     remap);
}

void PackageSet::JoinEdges(uint32_t pkg,
                           GraphTile& tile,
                           const GraphId& real_base,
                           NodeOwners& owners) const {
  const GraphTileHeader* header = tile.header_;
  const auto base_ll = header->base_ll();
  const int self = static_cast<int>(pkg);
  const auto disable = [](DirectedEdge& edge) {
    edge.set_forwardaccess(0);
    edge.set_reverseaccess(0);
    edge.shortcut_ = 0;
  };

  // Edges starting at nodes we do not own are never reached forward, and the edge bins drop them,
  // but the reverse of an edge crossing from one of our nodes can still be used: loki swaps in the
  // opposing edge of a disabled edge, and reverse searches only check the opposing edge's access.
  // Those reverse edges are disabled. With a single owner for the whole tile, an edge that stays
  // in it cannot end at a node of ours; with a foreign owner for everything the tile's edges reach
  // (from the reference index), no edge can, and the edges are not even read.
  const int tile_owner = BoxOwner(TileHierarchy::GetGraphIdBoundingBox(real_base));
  const bool check_foreign_edges = !ReachesOnlyForeign(pkg, real_base.value, tile_owner);

  // NodeLL and FindOwner of an edge's end node; this tile's nodes come from the owners found
  const auto end_node = [&](const GraphId& end_real, PointLL& end, int& owner) {
    if (end_real.tile_base() != real_base) {
      if (!NodeLL(pkg, end_real, end)) {
        return false;
      }
      owner = FindOwner(end);
      return true;
    }
    if (end_real.id() >= header->nodecount()) {
      return false;
    }
    end = tile.node(end_real.id())->latlng(base_ll);
    owner = owners.Of(end_real.id());
    return true;
  };

  uint64_t disabled = 0, joined = 0, lost = 0;
  for (uint32_t n = 0; n < header->nodecount(); ++n) {
    const NodeInfo* node = tile.node(n);
    const PointLL start = node->latlng(base_ll);
    if (owners.Of(n) != self) {
      disabled += node->edge_count();
      for (uint32_t k = 0; check_foreign_edges && k < node->edge_count(); ++k) {
        DirectedEdge& target = tile.directededges_[node->edge_index() + k];
        PointLL end;
        int end_owner = 0;
        if ((tile_owner != kUndecided && !target.leaves_tile()) ||
            !end_node(target.endnode(), end, end_owner) || end_owner != self) {
          continue;
        }
        DirectedEdge edge = target;
        edge.set_endnode(Virtual(pkg, edge.endnode()));
        disable(edge);
        Store(target, edge);
      }
      continue;
    }
    // not_thru was computed on this package's graph, which is truncated beyond its polygon; a
    // tile the package wholly owns keeps it
    const bool near_boundary = tile_owner != self && NearOwnershipBoundary(start);
    for (uint32_t k = 0; k < node->edge_count(); ++k) {
      DirectedEdge& target = tile.directededges_[node->edge_index() + k];
      DirectedEdge edge = target;
      const GraphId end_real = edge.endnode();
      edge.set_endnode(Virtual(pkg, end_real));
      if (near_boundary) {
        edge.set_not_thru(false);
      }
      PointLL end;
      int end_owner = 0;
      if (!end_node(end_real, end, end_owner)) {
        // the package does not have the end node: keep the edge out of routing
        if (RawTile(pkg, end_real.tile_base())) {
          LOG_WARN("Package {} edge {} ends at missing node {}", packages_[pkg]->name,
                   std::to_string(
                       GraphId(real_base.tileid(), real_base.level(), node->edge_index() + k)),
                   std::to_string(end_real));
        }
        edge.set_endnode(Nowhere(end_real));
        disable(edge);
        Store(target, edge);
        continue;
      }
      if (end_owner == self) {
        Store(target, edge);
        continue;
      }
      Twin twin{};
      if (!edge.is_shortcut() && FindTwin(pkg, tile, target, real_base, n, start, end,
                                          static_cast<uint32_t>(end_owner), twin)) {
        edge.set_endnode(twin.endnode);
        edge.set_opp_index(twin.opp_index);
        edge.set_opp_local_idx(twin.opp_local_idx);
        edge.set_restrictions(twin.restrictions);
        edge.set_reverseaccess(twin.reverseaccess);
        edge.set_not_thru(false);
        edge.set_leaves_tile(true);
        ++joined;
      } else {
        // crossing shortcuts are disabled too, with their shortcut mask
        disable(edge);
        if (!edge.is_shortcut()) {
          ++lost;
          LOG_DEBUG("Join lost: package {} level {} way {} length {} {},{} -> {},{} owned by {} ({})",
                    packages_[pkg]->name, real_base.level(), tile.edgeinfo(&target).wayid(),
                    edge.length(), start.lng(), start.lat(), end.lng(), end.lat(),
                    packages_[end_owner]->name, twin.reason);
        }
      }
      Store(target, edge);
    }
  }

  counters_.joined += joined;
  counters_.lost += lost;
  counters_.disabled += disabled;
}

bool PackageSet::ReachesOnlyForeign(uint32_t pkg, uint64_t tile_base, int tile_owner) const {
  if (tile_owner == kUndecided || tile_owner == static_cast<int>(pkg)) {
    return false;
  }
  const auto& package = *packages_[pkg];
  const auto refs = package.refs.find(tile_base);
  if (refs == package.refs.end()) {
    return false;
  }
  if (!refs->second.has_reach) {
    return true;
  }
  const int reach_owner = BoxOwner(refs->second.reach);
  return reach_owner != kUndecided && reach_owner != static_cast<int>(pkg);
}

bool PackageSet::SelfClean(uint32_t pkg, const TileEntry& entry) const {
  return entry.copy == 0 &&
         BoxOwner(TileHierarchy::GetGraphIdBoundingBox(GraphId(entry.real))) == static_cast<int>(pkg);
}

PackageSet::Work PackageSet::ComputeWork(uint32_t pkg, const TileEntry& entry) const {
  const auto& package = *packages_[pkg];
  const auto refs = package.refs.find(entry.real);
  if (refs == package.refs.end() || !SelfClean(pkg, entry)) {
    return Work::kFull;
  }
  // the nodes outside the tile that its edges, bins, and restrictions lead to must all be ours
  if (refs->second.has_reach && BoxOwner(refs->second.reach) != static_cast<int>(pkg)) {
    return Work::kFull;
  }
  // otherwise only ids change, and only if they point into a tile this package does not hold in
  // slot 0; references into tiles it does not ship at all need the full rewrite to disable them
  Work work = Work::kNone;
  const auto level = GraphId(entry.real).level();
  for (const auto ref : refs->second.ids) {
    const auto* other = package.Find(ref);
    if (!other) {
      return Work::kFull;
    }
    if (other->copy != 0) {
      if (GraphId(ref).level() == level) {
        return Work::kIdsOnly;
      }
      work = Work::kLinksOnly;
    }
  }
  return work;
}

graph_tile_ptr PackageSet::LoadTile(const GraphId& tile_base) const {
  const auto found = by_copy_.find(tile_base.tile_base().value);
  if (found == by_copy_.end()) {
    return nullptr;
  }
  const uint32_t pkg = found->second.first;
  const auto& package = *packages_[pkg];
  const TileEntry& entry = package.tiles[found->second.second];
  const GraphId real(entry.real);

  if (overlays_) {
    if (entry.overlay_size == 0) {
      if (!entry.counted.exchange(true)) {
        ++counters_.clean_tiles;
      }
      RawTile(pkg, real);
      return entry.raw;
    }
    std::call_once(entry.overlay_once, [&]() {
      if (entry.state == TileEntry::State::kBroken) {
        throw std::runtime_error("Tile " + std::to_string(real) + " of package " + package.name +
                                 " could not be patched");
      }
      const char* patch = package.overlay->data() + entry.overlay_offset;
      try {
        // the patch is applied inside the package's private mapping: only the pages it changes
        // become private memory
        std::call_once(package.private_once,
                       [&]() { package.private_map = std::make_shared<PrivateMap>(*package.file); });
        entry.overlay_tile =
            GraphTile::Create(GraphId(entry.tileid, real.level(), 0),
                              std::make_unique<MappedTileMemory<PrivateMap>>(package.private_map,
                                                                             package.private_map
                                                                                     ->data +
                                                                                 entry.offset,
                                                                             entry.size));
        // the layout is read from the tile as the package built it; the patch then changes it
        ApplyTilePatch(patch, entry.overlay_size, package.private_map->data + entry.offset,
                       entry.size, Layout(*entry.overlay_tile));
      } catch (...) {
        entry.state = TileEntry::State::kBroken;
        throw;
      }
      ++counters_.overlay_tiles;
    });
    return entry.overlay_tile;
  }

  // The runtime join: the tile is rewritten in place inside the package's private mapping, once.
  // Only this tile's lock is held while it is rewritten, so different tiles join in parallel; the
  // rewrite reads other tiles through the read-only mapping only.
  std::lock_guard<std::mutex> lock(entry.mutex);
  if (entry.state == TileEntry::State::kPending) {
    const Work work = ComputeWork(pkg, entry);
    if (work == Work::kNone) {
      entry.state = TileEntry::State::kClean;
      ++counters_.clean_tiles;
    } else {
      std::call_once(package.private_once,
                     [&]() { package.private_map = std::make_shared<PrivateMap>(*package.file); });
      auto tile =
          GraphTile::Create(real,
                            std::make_unique<MappedTileMemory<PrivateMap>>(package.private_map,
                                                                           package.private_map->data +
                                                                               entry.offset,
                                                                           entry.size));
      try {
        Rewrite(pkg, const_cast<GraphTile&>(*tile), work);
      } catch (...) {
        // part of the tile may be rewritten already: never serve it
        entry.state = TileEntry::State::kBroken;
        throw;
      }
      entry.joined = std::move(tile);
      entry.state = TileEntry::State::kRewritten;
      ++(work == Work::kFull ? counters_.full_rewrites : counters_.id_rewrites);
    }
  }
  switch (entry.state) {
    case TileEntry::State::kClean:
      RawTile(pkg, real);
      return entry.raw;
    case TileEntry::State::kRewritten:
      return entry.joined;
    default:
      throw std::runtime_error("Tile " + std::to_string(real) + " of package " + package.name +
                               " could not be joined");
  }
}

void PackageSet::ReleasePages() const {
#ifndef _WIN32
  for (const auto& package : packages_) {
    madvise(package->tar->mm.get(), package->tar->mm.size(), MADV_DONTNEED);
  }
#endif
}

PackageSet::JoinReport PackageSet::WriteOverlays(const boost::property_tree::ptree& pt,
                                                 const std::string& dir) {
  const auto packages = pt.get_child_optional("packages");
  if (!packages || packages->empty()) {
    throw std::runtime_error("No packages configured (mjolnir.packages)");
  }
  // not registered: the overlays get the copy ids of the set alone, whatever else the process holds
  auto set = Build(pt);
  set->IndexCopies();
  const auto ids = set->CopyIds();

  // mjolnir.concurrency, the thread count of Valhalla's other tools, defaults to every core
  const size_t threads = std::max<size_t>(
      1, pt.get<size_t>("concurrency", std::max(1u, std::thread::hardware_concurrency())));

  JoinReport report;
  report.key = set->join_key_;
  std::filesystem::create_directories(dir);
  std::vector<std::pair<std::string, std::string>> files; // (partial, final)
  try {
    // the patches of every package, in tile order, written once all are made
    std::vector<std::vector<std::pair<uint32_t, TilePatch>>> patches(set->packages_.size());
    std::vector<Overlay> overlays;
    for (uint32_t pkg = 0; pkg < set->packages_.size(); ++pkg) {
      const auto& package = *set->packages_[pkg];
      const auto file = (std::filesystem::path(dir) / (package.name + kOverlaySuffix)).string();
      files.emplace_back(file + ".partial", file);
      overlays.push_back(Overlay{package.name});
#if defined(POSIX_FADV_SEQUENTIAL)
      posix_fadvise(package.file->fd(), 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
    }

    // The tiles that need work, package by package in tile order. They are rewritten on
    // `threads` threads, each into a heap buffer read from its tar, and their patches are
    // collected here in this order, so the overlays don't depend on the thread count. At most a
    // window of tiles is in flight, and the packages' mapped pages are dropped every
    // kReleaseBytes: memory stays bounded by a few tiles per thread.
    struct Job {
      uint32_t pkg;
      size_t tile;
      Work work;
    };
    std::vector<Job> jobs;
    for (uint32_t pkg = 0; pkg < set->packages_.size(); ++pkg) {
      const auto& package = *set->packages_[pkg];
      for (size_t i = 0; i < package.tile_count; ++i) {
        const Work work = set->ComputeWork(pkg, package.tiles[i]);
        if (work == Work::kNone) {
          ++set->counters_.clean_tiles;
        } else {
          jobs.push_back({pkg, i, work});
        }
      }
    }
    // the patch of a rewritten tile (empty when the join left its bytes as they were), or why it
    // couldn't be made
    struct Done {
      bool ready = false;
      TilePatch patch;
      std::exception_ptr error;
    };
    const size_t window = 2 * threads;
    std::vector<Done> done(jobs.size());
    std::mutex mutex;
    std::condition_variable ready, room;
    size_t next = 0, written = 0;
    bool stop = false;
    const auto rewrite = [&]() {
      for (;;) {
        size_t k;
        {
          std::unique_lock<std::mutex> lock(mutex);
          room.wait(lock, [&] { return stop || next >= jobs.size() || next < written + window; });
          if (stop || next >= jobs.size()) {
            return;
          }
          k = next++;
        }
        Done result;
        try {
          const auto& job = jobs[k];
          const auto& package = *set->packages_[job.pkg];
          const auto& entry = package.tiles[job.tile];
          // pread, not the mapping: the tiles come in file order, so the reads are sequential,
          // while faulting the mapping in reads the file a few pages at a time
          std::vector<char> original(entry.size);
          package.file->Read(original.data(), entry.size, entry.offset);
          auto tile = GraphTile::Create(GraphId(entry.real), std::vector<char>(original));
          set->Rewrite(job.pkg, const_cast<GraphTile&>(*tile), job.work);
          ++(job.work == Work::kFull ? set->counters_.full_rewrites : set->counters_.id_rewrites);
          const auto* joined = reinterpret_cast<const char*>(tile->header());
          // only tiles whose bytes changed go into the overlay
          result.patch =
              MakeTilePatch(original.data(), joined, entry.size, Layout(*tile), kPatchLevel);
        } catch (...) {
          result.error = std::current_exception();
        }
        {
          std::lock_guard<std::mutex> lock(mutex);
          result.ready = true;
          done[k] = std::move(result);
        }
        ready.notify_all();
      }
    };
    std::vector<std::thread> workers;
    const auto join_workers = [&]() {
      {
        std::lock_guard<std::mutex> lock(mutex);
        stop = true;
      }
      room.notify_all();
      for (auto& worker : workers) {
        worker.join();
      }
      workers.clear();
    };
    try {
      for (size_t t = 0; t < std::min(threads, jobs.size()); ++t) {
        workers.emplace_back(rewrite);
      }
      size_t read = 0;
      for (size_t k = 0; k < jobs.size(); ++k) {
        Done result;
        {
          std::unique_lock<std::mutex> lock(mutex);
          ready.wait(lock, [&] { return done[k].ready; });
          result = std::move(done[k]);
          done[k] = Done{};
        }
        if (result.error) {
          std::rethrow_exception(result.error);
        }
        const auto& job = jobs[k];
        const auto& entry = set->packages_[job.pkg]->tiles[job.tile];
        if (!result.patch.bytes.empty()) {
          patches[job.pkg].emplace_back(static_cast<uint32_t>(entry.real), std::move(result.patch));
        }
        read += entry.size;
        if (read >= kReleaseBytes) {
          set->ReleasePages();
          read = 0;
        }
        {
          std::lock_guard<std::mutex> lock(mutex);
          written = k + 1;
        }
        room.notify_all();
      }
      join_workers();
    } catch (...) {
      join_workers();
      throw;
    }
    // durable before any is renamed into place; syncing them together lets the writes of one
    // package overlap the work on the next
    std::vector<std::unique_ptr<DurableFile>> outs;
    for (uint32_t pkg = 0; pkg < set->packages_.size(); ++pkg) {
      std::string directory, blobs, file;
      Hasher content;
      for (const auto& [tile, patch] : patches[pkg]) {
        PutU32(directory, tile);
        PutU32(directory, static_cast<uint32_t>(blobs.size()));
        PutU32(directory, static_cast<uint32_t>(patch.bytes.size()));
        blobs += patch.bytes;
        content.U64(tile);
        content.U64(patch.content_hash);
        if (blobs.size() > std::numeric_limits<uint32_t>::max()) {
          throw std::runtime_error("The overlay of " + set->packages_[pkg]->name + " is too large");
        }
      }
      for (const auto& c : ids) {
        PutU32(file, c.level);
        PutU32(file, c.tileid);
        PutU32(file, c.copy);
        PutU32(file, c.id);
      }
      file += directory;
      file += blobs;
      Hasher hasher;
      hasher.Bytes(file.data(), file.size());
      std::string header;
      PutU32(header, kOverlayMagic);
      PutU32(header, kOverlayVersion);
      PutU64(header, set->join_key_);
      PutU32(header, static_cast<uint32_t>(patches[pkg].size()));
      PutU32(header, static_cast<uint32_t>(ids.size()));
      PutU64(header, hasher.value());
      file.insert(0, header);
      outs.push_back(std::make_unique<DurableFile>(files[pkg].first));
      outs.back()->Write(file);
      overlays[pkg].tiles = patches[pkg].size();
      overlays[pkg].bytes = file.size();
      overlays[pkg].content_hash = content.value();
    }
    for (auto& out : outs) {
      out->Sync();
    }
    report.overlays = std::move(overlays);
    // all overlays are complete: a reader sees the old ones or the new ones, and a mix of both
    // has different join keys, which falls back to the runtime join
    for (const auto& [partial, file] : files) {
      std::filesystem::rename(partial, file);
    }
    SyncDirectory(dir);
  } catch (...) {
    for (const auto& file : files) {
      std::error_code ignored;
      std::filesystem::remove(file.first, ignored);
    }
    throw;
  }
  report.stats = set->stats();
  return report;
}

uint64_t PackageSet::WriteTileRefs(const std::string& tile_extract, const std::string& file) {
  auto archive = std::make_shared<midgard::tar>(tile_extract);
  std::vector<TarTile> list;
  const uint64_t fingerprint = ScanTar(*archive, tile_extract, list);
  std::unordered_map<uint64_t, graph_tile_ptr> tiles;
  for (const auto& t : list) {
    tiles.emplace(t.id,
                  GraphTile::Create(GraphId(t.id),
                                    std::make_unique<MappedTileMemory<midgard::tar>>(archive,
                                                                                     archive->mm
                                                                                             .get() +
                                                                                         t.offset,
                                                                                     t.size)));
  }
  const auto get_tile = [&](const GraphId& id) -> const GraphTile& {
    const auto tile = tiles.find(id.tile_base().value);
    if (tile == tiles.end()) {
      throw std::runtime_error(tile_extract + " references missing tile " +
                               std::to_string(id.tile_base()));
    }
    return *tile->second;
  };

  std::string out;
  PutU32(out, kTileRefsMagic);
  PutU32(out, kTileRefsVersion);
  PutU64(out, fingerprint);
  PutU32(out, static_cast<uint32_t>(list.size()));
  uint64_t total_refs = 0;
  // in id order, so the index is the same on every platform
  for (const auto& t : list) {
    const GraphId id(t.id);
    const auto& tile = get_tile(id);
    std::set<uint32_t> refs;
    float box[4] = {180.f, 90.f, -180.f, -90.f};
    const auto add = [&](const GraphId& ref) {
      if (ref.is_valid() && ref.tile_base() != id) {
        refs.insert(static_cast<uint32_t>(ref.tile_base().value));
      }
    };
    const auto reach = [&](const GraphId& node) {
      if (!node.is_valid() || node.tile_base() == id) {
        return;
      }
      const auto& other = get_tile(node);
      if (node.id() >= other.header()->nodecount()) {
        throw std::runtime_error(tile_extract + " references missing node " + std::to_string(node));
      }
      const auto ll = other.node(node.id())->latlng(other.header()->base_ll());
      box[0] = std::min(box[0], static_cast<float>(ll.lng()) - kReachPadDeg);
      box[1] = std::min(box[1], static_cast<float>(ll.lat()) - kReachPadDeg);
      box[2] = std::max(box[2], static_cast<float>(ll.lng()) + kReachPadDeg);
      box[3] = std::max(box[3], static_cast<float>(ll.lat()) + kReachPadDeg);
    };
    // an edge in another tile: its tile, and its start node
    const auto add_edge = [&](const GraphId& edge) {
      add(edge);
      if (edge.is_valid() && edge.tile_base() != id) {
        const auto& other = get_tile(edge);
        const int64_t start =
            edge.id() < other.header()->directededgecount() ? StartNode(other, edge.id()) : -1;
        if (start < 0) {
          throw std::runtime_error(
              tile_extract + " references an edge without a start node: " + std::to_string(edge));
        }
        reach(GraphId(edge.tileid(), edge.level(), static_cast<uint32_t>(start)));
      }
    };
    const auto* header = tile.header();
    for (uint32_t i = 0; i < header->directededgecount(); ++i) {
      const auto end = tile.directededge(i)->endnode();
      add(end);
      reach(end);
    }
    for (uint32_t i = 0; i < header->transitioncount(); ++i) {
      add(tile.transition(i)->endnode());
    }
    for (size_t b = 0; b < kBinCount; ++b) {
      for (const auto& edge : tile.GetBin(b)) {
        add_edge(edge);
      }
    }
    // every restriction LoadTile remaps, walked the same way
    const auto add_restriction = [&](ComplexRestriction& restriction) {
      add_edge(restriction.from_graphid());
      add_edge(restriction.to_graphid());
      restriction.WalkVias([&](const GraphId* via) {
        add_edge(*via);
        return WalkingVia::KeepWalking;
      });
    };
    ForEachRestriction(tile.complex_restriction_forward_, tile.complex_restriction_forward_size_,
                       add_restriction);
    ForEachRestriction(tile.complex_restriction_reverse_, tile.complex_restriction_reverse_size_,
                       add_restriction);

    PutU32(out, static_cast<uint32_t>(id.value));
    for (const float value : box) {
      uint32_t bits = 0;
      std::memcpy(&bits, &value, sizeof(bits));
      PutU32(out, bits);
    }
    PutU32(out, static_cast<uint32_t>(refs.size()));
    for (const auto ref : refs) {
      PutU32(out, ref);
    }
    total_refs += refs.size();
  }

  const auto partial = file + ".partial";
  {
    std::ofstream stream(partial, std::ios::binary | std::ios::trunc);
    stream.write(out.data(), static_cast<std::streamsize>(out.size()));
    stream.close();
    if (!stream) {
      std::remove(partial.c_str());
      throw std::runtime_error("Could not write " + file);
    }
  }
  std::filesystem::rename(partial, file);
  return total_refs;
}

} // namespace baldr
} // namespace valhalla
