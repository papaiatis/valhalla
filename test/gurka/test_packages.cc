// Joins independently built routing packages (mjolnir.packages) at runtime.
#include "baldr/graphreader.h"
#include "baldr/graphtileheader.h"
#include "baldr/packageset.h"
#include "baldr/tilehierarchy.h"
#include "gurka.h"
#include "microtar.h"
#include "midgard/sequence.h"
#include "test.h"

#include <boost/property_tree/ptree.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <map>
#include <random>
#include <set>
#include <string>
#include <thread>
#include <vector>

using namespace valhalla;
using namespace valhalla::baldr;
namespace fs = std::filesystem;

namespace {

constexpr double kJoinTolerance = 10;
constexpr double kMetersPerDegLat = 110574.0;
constexpr double kMetersPerDegLng = 111320.0;

std::string ReadFile(const std::string& file) {
  std::ifstream in(file, std::ios::binary);
  return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

void WriteFile(const std::string& file, const std::string& data) {
  std::ofstream(file, std::ios::binary | std::ios::trunc) << data;
}

// Writes the .gph tiles below tile_dir, except the excluded ones, into a tar, in path order. The
// creation date in the tile headers is cleared, so a fixture's tiles are the same bytes on every
// day and platform.
void TarTiles(const std::string& tile_dir,
              const std::string& tar_file,
              const std::set<GraphId>& excluded = {}) {
  std::vector<fs::path> tiles;
  for (const auto& entry : fs::recursive_directory_iterator(tile_dir)) {
    if (entry.is_regular_file() && entry.path().extension() == ".gph" &&
        !excluded.count(GraphId::FromTilePath(fs::relative(entry.path(), tile_dir).string()))) {
      tiles.push_back(entry.path());
    }
  }
  std::sort(tiles.begin(), tiles.end());
  mtar_t tar;
  ASSERT_EQ(mtar_open(&tar, tar_file.c_str(), "w"), MTAR_ESUCCESS);
  for (const auto& tile : tiles) {
    auto data = ReadFile(tile.string());
    ASSERT_GE(data.size(), sizeof(GraphTileHeader));
    GraphTileHeader header;
    std::memcpy(&header, data.data(), sizeof(header));
    header.set_date_created(0);
    std::memcpy(data.data(), &header, sizeof(header));
    const auto name = fs::relative(tile, tile_dir).generic_string();
    ASSERT_EQ(mtar_write_file_header(&tar, name.c_str(), static_cast<unsigned>(data.size())),
              MTAR_ESUCCESS);
    ASSERT_EQ(mtar_write_data(&tar, data.data(), static_cast<unsigned>(data.size())), MTAR_ESUCCESS);
  }
  mtar_finalize(&tar);
  mtar_close(&tar);
}

void WriteBox(const std::string& file,
              double min_lng,
              double min_lat,
              double max_lng,
              double max_lat) {
  std::ofstream out(file);
  out << "region\nring\n"
      << min_lng << " " << min_lat << "\n"
      << max_lng << " " << min_lat << "\n"
      << max_lng << " " << max_lat << "\n"
      << min_lng << " " << max_lat << "\n"
      << "END\nEND\n";
}

// A package built for a test: its files and the build time of its mjolnir.packages entry.
struct Package {
  std::string name, tar, poly, refs;
  int64_t build_time = 0;
};

// The ways with a node whose longitude lies in [min_lng, max_lng].
gurka::ways
WaysWithin(const gurka::ways& ways, const gurka::nodelayout& layout, double min_lng, double max_lng) {
  gurka::ways within;
  for (const auto& way : ways) {
    if (std::any_of(way.first.begin(), way.first.end(), [&](char node) {
          const double lng = layout.at(std::string(1, node)).lng();
          return lng >= min_lng && lng <= max_lng;
        })) {
      within.insert(way);
    }
  }
  return within;
}

// The relations all of whose ways are among the given ones, as an extract would keep them.
gurka::relations RelationsWithin(const gurka::relations& relations, const gurka::ways& ways) {
  gurka::relations within;
  for (const auto& relation : relations) {
    if (std::all_of(relation.members.begin(), relation.members.end(), [&](const auto& member) {
          return member.type != gurka::way_member || ways.count(member.ref);
        })) {
      within.push_back(relation);
    }
  }
  return within;
}

// Pins an OSM id on every way, so each package numbers the same way the same.
gurka::ways WithWayIds(gurka::ways ways) {
  uint64_t id = 1000;
  for (auto& way : ways) {
    way.second.emplace("osm_id", std::to_string(id++));
  }
  return ways;
}

// Builds a package's tiles from the given ways; its polygon is the box given.
// Tiles of a package's build left out of its tar leave its references into them dangling, which
// also keeps it from having a tile reference index.
Package BuildPackage(const std::string& dir,
                     const std::string& name,
                     const gurka::nodelayout& layout,
                     const gurka::ways& ways,
                     const gurka::relations& relations,
                     const midgard::AABB2<midgard::PointLL>& polygon,
                     int64_t build_time = 0,
                     const std::set<GraphId>& excluded = {}) {
  const auto workdir = dir + "/" + name;
  gurka::buildtiles(layout, ways, {}, RelationsWithin(relations, ways), workdir);
  Package package{name, dir + "/" + name + ".tar", dir + "/" + name + ".poly",
                  excluded.empty() ? dir + "/" + name + "_tile_refs.bin" : "", build_time};
  TarTiles(workdir, package.tar, excluded);
  WriteBox(package.poly, polygon.minx(), polygon.miny(), polygon.maxx(), polygon.maxy());
  if (!package.refs.empty()) {
    PackageSet::WriteTileRefs(package.tar, package.refs);
  }
  return package;
}

// The mjolnir.packages list of the packages.
boost::property_tree::ptree PackageList(const std::vector<Package>& packages) {
  boost::property_tree::ptree list;
  for (const auto& package : packages) {
    boost::property_tree::ptree entry;
    entry.put("name", package.name);
    entry.put("tile_extract", package.tar);
    entry.put("polygon", package.poly);
    entry.put("build_time", package.build_time);
    if (!package.refs.empty()) {
      entry.put("tile_refs", package.refs);
    }
    list.push_back({"", entry});
  }
  return list;
}

// A routable map over the packages.
gurka::map PackageMap(const std::vector<Package>& packages,
                      const std::string& dir,
                      const gurka::nodelayout& layout,
                      const std::string& overlays = "") {
  gurka::map map{test::make_config(dir + "/config"), layout};
  // tiles come from the packages only
  for (const auto* key : {"tile_extract", "traffic_extract"}) {
    map.config.get_child("mjolnir").erase(key);
  }
  map.config.put("mjolnir.tile_dir", dir + "/no_tiles");
  map.config.put_child("mjolnir.packages", PackageList(packages));
  map.config.put("mjolnir.package_join_tolerance", kJoinTolerance);
  if (!overlays.empty()) {
    map.config.put("mjolnir.package_joined", overlays);
  }
  return map;
}

// The same map, routed through the overlays in the given directory.
gurka::map WithOverlays(gurka::map map, const std::string& overlays) {
  map.config.put("mjolnir.package_joined", overlays);
  return map;
}

// A route's edges as (way id, length) and its totals: equal for equal routes over graphs that
// number their edges differently.
struct RouteSummary {
  std::vector<std::pair<uint64_t, uint32_t>> edges;
  double km = 0, seconds = 0;
};

RouteSummary Route(const gurka::map& map, const std::string& from, const std::string& to) {
  const auto result = gurka::do_action(Options::route, map, {from, to}, "auto");
  RouteSummary route;
  const auto& leg = result.trip().routes(0).legs(0);
  for (int i = 0; i < leg.node_size(); ++i) {
    if (leg.node(i).has_edge()) {
      route.edges.emplace_back(leg.node(i).edge().way_id(),
                               static_cast<uint32_t>(
                                   std::lround(leg.node(i).edge().length_km() * 1000)));
    }
  }
  route.km = result.directions().routes(0).legs(0).summary().length();
  route.seconds = result.directions().routes(0).legs(0).summary().time();
  return route;
}

// The way ids along a route.
std::vector<uint64_t>
RouteWays(const gurka::map& map, const std::string& from, const std::string& to) {
  std::vector<uint64_t> ways;
  for (const auto& edge : Route(map, from, to).edges) {
    if (ways.empty() || ways.back() != edge.first) {
      ways.push_back(edge.first);
    }
  }
  return ways;
}

uint64_t WayId(const gurka::ways& ways, const std::string& name) {
  return std::stoull(ways.at(name).at("osm_id"));
}

// Straight-line distance in km between two nodes of the layout.
double Km(const gurka::nodelayout& layout, const std::string& from, const std::string& to) {
  return layout.at(from).Distance(layout.at(to)) / 1000;
}

// Every tile's bytes, by tile id.
std::map<GraphId, std::string> TileBytes(GraphReader& reader) {
  std::map<GraphId, std::string> bytes;
  for (const auto& id : reader.GetTileSet()) {
    const auto tile = reader.GetGraphTile(id);
    if (tile) {
      bytes[id] =
          std::string(reinterpret_cast<const char*>(tile->header()), tile->header()->end_offset());
    }
  }
  return bytes;
}

// Every usable edge whose start node its package owns must have an opposing edge leading back to
// its start node, across package joins too: a reference into another package's numbering breaks
// this. Edges starting at nodes another package owns are never reached, and are not rewritten.
void ExpectConsistentGraph(GraphReader& reader) {
  const auto& set = *reader.package_set();
  size_t checked = 0;
  for (const auto& tile_id : reader.GetTileSet()) {
    const auto tile = reader.GetGraphTile(tile_id);
    ASSERT_TRUE(tile);
    const int package = set.PackageOf(tile_id);
    ASSERT_GE(package, 0);
    for (uint32_t n = 0; n < tile->header()->nodecount(); ++n) {
      const auto* node = tile->node(n);
      if (set.Owner(node->latlng(tile->header()->base_ll())) != package) {
        continue;
      }
      for (uint32_t k = 0; k < node->edge_count(); ++k) {
        const GraphId edge_id(tile_id.tileid(), tile_id.level(), node->edge_index() + k);
        const auto* edge = tile->directededge(edge_id);
        if (!(edge->forwardaccess() & kAutoAccess)) {
          continue;
        }
        graph_tile_ptr opp_tile;
        const DirectedEdge* opp = nullptr;
        try {
          opp = reader.GetOpposingEdge(edge_id, opp_tile);
        } catch (const std::exception& e) {
          ADD_FAILURE() << "opposing edge of " << edge_id << ": " << e.what();
          continue;
        }
        ASSERT_NE(opp, nullptr) << "no opposing edge of " << edge_id;
        EXPECT_EQ(opp->endnode(), GraphId(tile_id.tileid(), tile_id.level(), n))
            << "opposing edge of " << edge_id << " does not lead back to its start node";
        ++checked;
      }
    }
  }
  EXPECT_GT(checked, 0u);
}

// The tiles of a tar, by tile id; the entries that are no tiles are left out.
std::map<GraphId, std::string> TarTileBytes(const std::string& file) {
  std::map<GraphId, std::string> tiles;
  midgard::tar archive(file);
  archive.for_each([&](const std::string& name, const char* data, size_t size) {
    try {
      tiles[GraphId::FromTilePath(name)] = std::string(data, size);
    } catch (const std::exception&) {}
    return true;
  });
  return tiles;
}

// Every tile of every package copy, by copy id, loaded through a package set.
std::map<GraphId, std::string> LoadAll(const PackageSet& set) {
  std::map<GraphId, std::string> tiles;
  for (const auto& id : set.AllTiles()) {
    const auto tile = set.LoadTile(id);
    EXPECT_TRUE(tile) << id;
    if (tile) {
      tiles[id] =
          std::string(reinterpret_cast<const char*>(tile->header()), tile->header()->end_offset());
    }
  }
  return tiles;
}

// The copy id of a package's copy of a real tile, or an invalid id.
GraphId CopyOf(const PackageSet& set, const std::string& package, const GraphId& real) {
  for (const auto& copy : set.Copies(real)) {
    if (set.name(set.PackageOf(copy)) == package) {
      return copy;
    }
  }
  return {};
}

// FNV-1a over a file's bytes: the same value on every platform.
uint64_t Fnv1a(const std::string& data) {
  uint64_t hash = 14695981039346656037ull;
  for (const unsigned char c : data) {
    hash = (hash ^ c) * 1099511628211ull;
  }
  return hash;
}

std::string OverlayFile(const std::string& dir, const std::string& package) {
  return dir + "/" + package + ".joined";
}

// ---------------------------------------------------------------------------------------------
// One road across a border between B and C. Package "west" is built from A-B-C, package "east"
// from B-C-D: both hold the way B-C, and every node lies in one level-2 tile.

const std::string kAscii = R"(
    A-------B-------C-------D
)";
constexpr double kGridSize = 250;
const midgard::PointLL kOrigin{5.1, 52.1};

// Builds both packages below dir; east_bc_way_id is the OSM id east uses for the way B-C.
std::vector<Package> BuildPackages(const std::string& dir, const std::string& east_bc_way_id) {
  const auto layout = gurka::detail::map_to_coordinates(kAscii, kGridSize, kOrigin);
  const double border = (layout.at("B").lng() + layout.at("C").lng()) / 2;
  const gurka::ways west_ways = {
      {"AB", {{"highway", "residential"}, {"osm_id", "1"}}},
      {"BC", {{"highway", "residential"}, {"osm_id", "2"}}},
  };
  const gurka::ways east_ways = {
      {"BC", {{"highway", "residential"}, {"osm_id", east_bc_way_id}}},
      {"CD", {{"highway", "residential"}, {"osm_id", "3"}}},
  };
  return {BuildPackage(dir, "west", layout, west_ways, {}, {4.9, 51.9, border, 52.3}),
          BuildPackage(dir, "east", layout, east_ways, {}, {border, 51.9, 5.4, 52.3})};
}

gurka::map JoinedMap(const std::vector<Package>& packages, const std::string& dir) {
  return PackageMap(packages, dir, gurka::detail::map_to_coordinates(kAscii, kGridSize, kOrigin));
}

double RouteKm(const gurka::map& map, const std::string& from, const std::string& to) {
  return Route(map, from, to).km;
}

double Km(const std::string& from, const std::string& to) {
  return Km(gurka::detail::map_to_coordinates(kAscii, kGridSize, kOrigin), from, to);
}

class PackageJoin : public ::testing::Test {
protected:
  static void SetUpTestSuite() {
    dir_ = VALHALLA_BUILD_DIR "test/data/gurka_packages/joined";
    fs::remove_all(dir_);
    fs::create_directories(dir_);
    packages_ = BuildPackages(dir_, "2");
    map_ = JoinedMap(packages_, dir_);
  }
  static std::string dir_;
  static std::vector<Package> packages_;
  static gurka::map map_;
};
std::string PackageJoin::dir_;
std::vector<Package> PackageJoin::packages_;
gurka::map PackageJoin::map_;

TEST_F(PackageJoin, RoutesAcrossTheJoin) {
  EXPECT_NEAR(RouteKm(map_, "A", "D"), Km("A", "D"), 0.01);
  EXPECT_NEAR(RouteKm(map_, "D", "A"), Km("D", "A"), 0.01);
  // B-C exists in both packages but is routed once
  EXPECT_NEAR(RouteKm(map_, "B", "C"), Km("B", "C"), 0.01);
}

TEST_F(PackageJoin, JoinCounts) {
  const auto set = PackageSet::FromConfig(map_.config.get_child("mjolnir"));
  ASSERT_TRUE(set);
  for (const auto& id : set->AllTiles()) {
    set->LoadTile(id);
  }
  const auto stats = set->stats();
  // B->C in west and C->B in east
  EXPECT_EQ(stats.joined, 2u);
  EXPECT_EQ(stats.lost, 0u);
  // west's C->B and east's B->C start at nodes of the other package
  EXPECT_EQ(stats.disabled, 2u);
}

TEST_F(PackageJoin, CopyIds) {
  GraphReader reader(map_.config.get_child("mjolnir"));
  const auto layout = gurka::detail::map_to_coordinates(kAscii, kGridSize, kOrigin);
  const auto real = TileHierarchy::GetGraphId(layout.at("B"), 2);
  const auto copies = reader.GetTileCopies(real);
  ASSERT_EQ(copies.size(), 2u);
  EXPECT_EQ(copies[0]->id(), real) << "slot 0 keeps the real id";
  const GraphId copy = copies[1]->id();
  EXPECT_NE(copy, real);
  EXPECT_GE(copy.tileid(), TileHierarchy::levels()[2].tiles.TileCount()) << "above the grid";
  EXPECT_EQ(PackageSet::Real(copy), real);
  EXPECT_EQ(PackageSet::CopyIndex(copy), 1u);
  EXPECT_EQ(PackageSet::Real(GraphId(copy.tileid(), 2, 7)), GraphId(real.tileid(), 2, 7));
  EXPECT_TRUE(copies[1]->BoundingBox() == copies[0]->BoundingBox());
  EXPECT_TRUE(reader.DoesTileExist(copy));
  EXPECT_EQ(reader.GetTileSet(2).size(), 2u);
  // the same copy of a tile keeps its id in every package set of the process
  const auto again = PackageSet::FromConfig(map_.config.get_child("mjolnir"));
  const auto ids = again->Copies(real);
  ASSERT_EQ(ids.size(), 2u);
  EXPECT_EQ(ids[1], copy);
}

TEST_F(PackageJoin, ConsistentGraph) {
  GraphReader reader(map_.config.get_child("mjolnir"));
  ExpectConsistentGraph(reader);
}

TEST_F(PackageJoin, SinglePackage) {
  auto map = map_;
  auto& list = map.config.get_child("mjolnir.packages");
  list.erase(list.begin());
  EXPECT_NEAR(RouteKm(map, "D", "C"), Km("D", "C"), 0.01);
  // east owns B beyond its polygon too, as the nearest package: A snaps to the end of B-C
  EXPECT_NEAR(RouteKm(map, "D", "A"), Km("D", "B"), 0.01);
}

TEST_F(PackageJoin, OverlaysRouteAndRewriteNothing) {
  const auto dir = dir_ + "/overlays";
  const auto report = PackageSet::WriteOverlays(map_.config.get_child("mjolnir"), dir);
  EXPECT_EQ(report.stats.joined, 2u);
  const auto map = WithOverlays(map_, dir);
  EXPECT_NEAR(RouteKm(map, "A", "D"), Km("A", "D"), 0.01);
  const auto set = PackageSet::FromConfig(map.config.get_child("mjolnir"));
  ASSERT_TRUE(set->uses_overlays());
  LoadAll(*set);
  EXPECT_EQ(set->stats().full_rewrites, 0u);
  EXPECT_EQ(set->stats().id_rewrites, 0u);
  GraphReader reader(map.config.get_child("mjolnir"));
  ExpectConsistentGraph(reader);
}

// Lines of this process's memory map that name the file.
size_t Mappings(const std::string& file) {
  size_t count = 0;
#ifdef __linux__
  std::ifstream maps("/proc/self/maps");
  const auto path = fs::canonical(file).string();
  for (std::string line; std::getline(maps, line);) {
    count += line.find(path) != std::string::npos;
  }
#endif
  return count;
}

TEST_F(PackageJoin, SetIsReleasedWithItsLastReader) {
  const auto config = map_.config.get_child("mjolnir");
  std::weak_ptr<const PackageSet> released;
  {
    GraphReader first(config);
    GraphReader second(config);
    ASSERT_TRUE(first.package_set());
    EXPECT_EQ(first.package_set(), second.package_set());
    released = first.package_set();
    for (const auto& id : first.GetTileSet()) {
      first.GetGraphTile(id);
    }
#ifdef __linux__
    EXPECT_GT(Mappings(packages_[0].tar), 0u);
#endif
  }
  EXPECT_TRUE(released.expired());
  EXPECT_EQ(Mappings(packages_[0].tar), 0u);
  GraphReader third(config);
  ASSERT_TRUE(third.package_set());
  EXPECT_NEAR(RouteKm(map_, "A", "D"), Km("A", "D"), 0.01);
}

TEST_F(PackageJoin, InvalidConfigsThrow) {
  auto config = map_.config.get_child("mjolnir");
  auto duplicate = config;
  duplicate.get_child("packages").back().second.put("name", "west");
  EXPECT_THROW(PackageSet::FromConfig(duplicate), std::runtime_error);

  for (const auto* name : {"", ".", "..", "a/b"}) {
    auto bad_name = config;
    bad_name.get_child("packages").back().second.put("name", name);
    EXPECT_THROW(PackageSet::FromConfig(bad_name), std::runtime_error) << name;
  }

  auto missing_tar = config;
  missing_tar.get_child("packages").back().second.put("tile_extract", dir_ + "/nonexistent.tar");
  EXPECT_THROW(PackageSet::FromConfig(missing_tar), std::runtime_error);

  auto missing_poly = config;
  missing_poly.get_child("packages").back().second.put("polygon", dir_ + "/nonexistent.poly");
  EXPECT_THROW(PackageSet::FromConfig(missing_poly), std::runtime_error);

  auto negative = config;
  negative.put("package_join_tolerance", -1);
  EXPECT_THROW(PackageSet::FromConfig(negative), std::runtime_error);

  // cut inside the last tile's data, at a tar block boundary
  const auto tar = ReadFile(packages_[0].tar);
  ASSERT_GT(tar.size(), 4 * 512u);
  WriteFile(dir_ + "/truncated.tar", tar.substr(0, tar.size() - 3 * 512));
  auto truncated = config;
  truncated.get_child("packages").front().second.put("tile_extract", dir_ + "/truncated.tar");
  EXPECT_THROW(PackageSet::FromConfig(truncated), std::runtime_error);
}

// A package tar replaced (renamed over) while a set uses it: the set keeps joining the file it
// loaded, whose tiles it indexed.
TEST_F(PackageJoin, ReplacedTarKeepsServingTheLoadedFile) {
  const auto dir = dir_ + "/replaced";
  fs::remove_all(dir);
  fs::create_directories(dir);
  auto config = map_.config.get_child("mjolnir");
  for (auto& entry : config.get_child("packages")) {
    const auto tar = dir + "/" + entry.second.get<std::string>("name") + ".tar";
    fs::copy_file(entry.second.get<std::string>("tile_extract"), tar);
    entry.second.put("tile_extract", tar);
  }
  const auto expected = LoadAll(*PackageSet::FromConfig(config));
  ASSERT_FALSE(expected.empty());

  const auto set = PackageSet::FromConfig(config);
  ASSERT_FALSE(set->uses_overlays());
  // another tar under both names: west's is east's, east's a copy with other tile data
  const auto west = dir + "/west.tar", east = dir + "/east.tar";
  auto other = ReadFile(east);
  for (size_t i = 3 * 512; i < other.size() - 2 * 512; ++i) {
    other[i] = static_cast<char>(~other[i]);
  }
  fs::copy_file(east, dir + "/next_west.tar");
  fs::rename(dir + "/next_west.tar", west);
  WriteFile(dir + "/next_east.tar", other);
  fs::rename(dir + "/next_east.tar", east);
  EXPECT_TRUE(LoadAll(*set) == expected);
}

// The join key follows the tile data even when tile headers hold no hash of it.
TEST_F(PackageJoin, JoinKeyHashesTilesWithoutDataHash) {
  const auto dir = dir_ + "/unhashed";
  fs::remove_all(dir);
  fs::create_directories(dir);
  // the west tar with the data hash cleared in every tile header
  auto tar = ReadFile(packages_[0].tar);
  size_t last_tile = 0;
  midgard::tar archive(packages_[0].tar);
  archive.for_each([&](const std::string&, const char* data, size_t size) {
    const size_t offset = static_cast<size_t>(data - archive.mm.get());
    GraphTileHeader header;
    std::memcpy(&header, tar.data() + offset, sizeof(header));
    header.set_raw_checksum(0);
    std::memcpy(&tar[offset], &header, sizeof(header));
    last_tile = offset + size - 1;
    return true;
  });
  auto config = map_.config.get_child("mjolnir");
  const auto key_of = [&](const std::string& data) {
    const auto file = dir + "/west.tar";
    WriteFile(file, data);
    config.get_child("packages").front().second.put("tile_extract", file);
    config.get_child("packages").front().second.erase("tile_refs");
    return PackageSet::FromConfig(config)->join_key();
  };
  const auto key = key_of(tar);
  EXPECT_EQ(key_of(tar), key);
  // a changed byte of tile data, header and size unchanged
  tar[last_tile] = static_cast<char>(~tar[last_tile]);
  EXPECT_NE(key_of(tar), key);
}

// east numbers the way B-C differently, so the packages cannot be joined there
class LostJoin : public ::testing::Test {
protected:
  static void SetUpTestSuite() {
    dir_ = VALHALLA_BUILD_DIR "test/data/gurka_packages/lost";
    fs::remove_all(dir_);
    fs::create_directories(dir_);
    map_ = JoinedMap(BuildPackages(dir_, "20"), dir_);
  }
  static std::string dir_;
  static gurka::map map_;
};
std::string LostJoin::dir_;
gurka::map LostJoin::map_;

TEST_F(LostJoin, NoRouteAcrossAndConsistentGraph) {
  const auto set = PackageSet::FromConfig(map_.config.get_child("mjolnir"));
  for (const auto& id : set->AllTiles()) {
    set->LoadTile(id);
  }
  EXPECT_EQ(set->stats().joined, 0u);
  EXPECT_EQ(set->stats().lost, 2u);

  EXPECT_THROW(RouteKm(map_, "A", "D"), std::exception);
  EXPECT_THROW(RouteKm(map_, "D", "A"), std::exception);
  EXPECT_NEAR(RouteKm(map_, "A", "B"), Km("A", "B"), 0.01);
  EXPECT_NEAR(RouteKm(map_, "D", "C"), Km("D", "C"), 0.01);
  // the reverse of a lost edge starts at a node its package does not own: it must be disabled,
  // or its unrewritten end node reads another package's tile
  GraphReader reader(map_.config.get_child("mjolnir"));
  ExpectConsistentGraph(reader);
}

// ---------------------------------------------------------------------------------------------
// One graph split into two packages, "west" and "east", by a line of longitude. Each package is
// built from the ways with a node within `reach` of its side, the way a region's extract reaches
// past its polygon: far enough that every edge attribute computed from the neighbourhood (road
// density within 2 km) is the same as in the single graph. The polygons overlap around the line
// like buffered region polygons do; with equal build times the line splits the overlap.

struct SplitMaps {
  gurka::map single;
  gurka::map joined;
  std::vector<Package> packages;
  gurka::ways ways;
};

SplitMaps BuildSplit(const std::string& dir,
                     const gurka::nodelayout& layout,
                     const gurka::ways& ways,
                     const gurka::relations& relations,
                     double line_lng,
                     double reach_m,
                     const gurka::nodelayout& east_layout = {}) {
  fs::remove_all(dir);
  fs::create_directories(dir);
  SplitMaps maps;
  maps.ways = WithWayIds(ways);
  maps.single = gurka::buildtiles(layout, maps.ways, {}, relations, dir + "/single");
  const double lat = layout.begin()->second.lat();
  const double reach = reach_m / (kMetersPerDegLng * std::cos(lat * M_PI / 180));
  const double overlap = reach / 10;
  const auto& east_nodes = east_layout.empty() ? layout : east_layout;
  maps.packages = {
      BuildPackage(dir, "west", layout, WaysWithin(maps.ways, layout, -180, line_lng + reach),
                   relations, {line_lng - 1, lat - 1, line_lng + overlap, lat + 1}),
      BuildPackage(dir, "east", east_nodes, WaysWithin(maps.ways, layout, line_lng - reach, 180),
                   relations, {line_lng - overlap, lat - 1, line_lng + 1, lat + 1}),
  };
  maps.joined = PackageMap(maps.packages, dir, layout);
  return maps;
}

// Rows of roads across the line, joined by north-south roads. The line runs between columns D and
// E. Row 1 is a secondary road (level 1), row 2 a primary road (level 0), the others residential
// (level 2); row 3 is two ways that meet next to the line. The north-south roads farther than
// `reach` from the line are in one package only, so the packages number their edges differently.
const std::vector<int> kSplitColumns = {0, 8, 16, 22, 28, 34, 42, 50};
const std::vector<std::string> kSplitRows = {"ABCDEFGH", "IJKLMNOP", "QRSTUVWX", "abcdefgh"};

std::string SplitAscii() {
  std::string ascii = "\n";
  for (size_t r = 0; r < kSplitRows.size(); ++r) {
    std::string row(kSplitColumns.back() + 1, '-'), gap(kSplitColumns.back() + 1, ' ');
    for (size_t c = 0; c < kSplitColumns.size(); ++c) {
      row[kSplitColumns[c]] = kSplitRows[r][c];
      gap[kSplitColumns[c]] = '|';
    }
    ascii += row + "\n" + (r + 1 < kSplitRows.size() ? gap + "\n" : "");
  }
  return ascii;
}

gurka::ways SplitWays() {
  gurka::ways ways = {
      {"ABCDEFGH", {{"highway", "residential"}}}, {"IJKLMNOP", {{"highway", "secondary"}}},
      {"QRSTUVWX", {{"highway", "primary"}}},     {"abcd", {{"highway", "residential"}}},
      {"defgh", {{"highway", "residential"}}},
  };
  for (size_t c = 0; c < kSplitColumns.size(); ++c) {
    std::string column;
    for (const auto& row : kSplitRows) {
      column += row[c];
    }
    ways[column] = {{"highway", "residential"}};
  }
  return ways;
}

class SplitJoin : public ::testing::Test {
protected:
  static void SetUpTestSuite() {
    // cells of about 120 m; the grid spans two level-2 tiles
    layout_ = gurka::detail::map_to_coordinates(SplitAscii(), 200, {5.2, 52.1});
    // a road in a level-2 tile of its own on each side, far from the line: tiles that need no
    // rewrite with a tile reference index
    auto ways = SplitWays();
    layout_["1"] = {4.9, 52.1};
    layout_["2"] = {4.91, 52.1};
    layout_["3"] = {5.65, 52.1};
    layout_["4"] = {5.66, 52.1};
    ways["12"] = {{"highway", "residential"}};
    ways["34"] = {{"highway", "residential"}};
    const double line = (layout_.at("D").lng() + layout_.at("E").lng()) / 2;
    maps_ = BuildSplit(VALHALLA_BUILD_DIR "test/data/gurka_packages/split", layout_, ways, {}, line,
                       3000);
  }
  static gurka::nodelayout layout_;
  static SplitMaps maps_;
};
gurka::nodelayout SplitJoin::layout_;
SplitMaps SplitJoin::maps_;

TEST_F(SplitJoin, NoLostJoins) {
  const auto set = PackageSet::FromConfig(maps_.joined.config.get_child("mjolnir"));
  for (const auto& id : set->AllTiles()) {
    set->LoadTile(id);
  }
  // one crossing edge per row and direction, on every level the rows are on
  EXPECT_EQ(set->stats().lost, 0u);
  EXPECT_GE(set->stats().joined, 2 * kSplitRows.size());
}

TEST_F(SplitJoin, RoutesEqualTheSingleGraph) {
  const std::vector<std::string> west = {"A", "J", "S", "d", "B", "L"};
  const std::vector<std::string> east = {"E", "O", "W", "h"};
  size_t routes = 0;
  for (const auto& w : west) {
    for (const auto& e : east) {
      for (const auto& [from, to] : {std::make_pair(w, e), std::make_pair(e, w)}) {
        const auto expected = Route(maps_.single, from, to);
        const auto actual = Route(maps_.joined, from, to);
        EXPECT_EQ(actual.edges, expected.edges) << from << " -> " << to;
        EXPECT_NEAR(actual.km, expected.km, 1e-3) << from << " -> " << to;
        EXPECT_NEAR(actual.seconds, expected.seconds, 1e-3) << from << " -> " << to;
        ++routes;
      }
    }
  }
  EXPECT_GE(routes, 20u);
}

TEST_F(SplitJoin, ConsistentGraph) {
  GraphReader reader(maps_.joined.config.get_child("mjolnir"));
  ExpectConsistentGraph(reader);
}

// Loading every tile from many threads at once gives the same tiles as loading them one by one.
// Readers on 8 threads, each loading every tile in its own order from one shared package set, get
// the tiles of a serial load.
void ExpectParallelLoadsEqual(const boost::property_tree::ptree& config,
                              const std::map<GraphId, std::string>& expected) {
  // every round's readers share one package set, released with them
  for (int round = 0; round < 3; ++round) {
    std::vector<std::map<GraphId, std::string>> loaded(8);
    std::vector<std::thread> threads;
    for (size_t t = 0; t < loaded.size(); ++t) {
      threads.emplace_back([&, t]() {
        GraphReader reader(config);
        auto ids = reader.GetTileSet();
        std::vector<GraphId> order(ids.begin(), ids.end());
        std::sort(order.begin(), order.end());
        std::shuffle(order.begin(), order.end(), std::mt19937(static_cast<unsigned>(t + round)));
        for (const auto& id : order) {
          const auto tile = reader.GetGraphTile(id);
          ASSERT_TRUE(tile);
          loaded[t][id] = std::string(reinterpret_cast<const char*>(tile->header()),
                                      tile->header()->end_offset());
        }
      });
    }
    for (auto& thread : threads) {
      thread.join();
    }
    for (const auto& tiles : loaded) {
      EXPECT_TRUE(tiles == expected);
    }
  }
}

TEST_F(SplitJoin, ParallelLoadsEqualSerialLoads) {
  const auto config = maps_.joined.config.get_child("mjolnir");
  // a package set of its own, not shared with the readers below
  const auto expected = LoadAll(*PackageSet::FromConfig(config));
  ASSERT_FALSE(expected.empty());
  ExpectParallelLoadsEqual(config, expected);
}

// The tiles of the overlays are patched in the packages' private mappings by whichever thread asks
// first: threads that share a page of a mapping must not disturb each other.
TEST_F(SplitJoin, ParallelLoadsFromOverlaysEqualSerialLoads) {
  const auto config = maps_.joined.config.get_child("mjolnir");
  const auto expected = LoadAll(*PackageSet::FromConfig(config));
  const std::string dir = VALHALLA_BUILD_DIR "test/data/gurka_packages/split/parallel_overlays";
  fs::remove_all(dir);
  PackageSet::WriteOverlays(config, dir);
  const auto with_overlays = WithOverlays(maps_.joined, dir).config.get_child("mjolnir");
  ASSERT_TRUE(PackageSet::FromConfig(with_overlays)->uses_overlays());
  ExpectParallelLoadsEqual(with_overlays, expected);
}

// The routes of RoutesEqualTheSingleGraph through a few of the pairs, for maps that must route
// like the single graph.
void ExpectSingleGraphRoutes(const gurka::map& map, const gurka::map& single) {
  for (const auto& [from, to] : {std::make_pair("A", "E"), std::make_pair("h", "S"),
                                 std::make_pair("J", "W"), std::make_pair("O", "d")}) {
    const auto expected = Route(single, from, to);
    const auto actual = Route(map, from, to);
    EXPECT_EQ(actual.edges, expected.edges) << from << " -> " << to;
    EXPECT_NEAR(actual.seconds, expected.seconds, 1e-3) << from << " -> " << to;
  }
}

TEST_F(SplitJoin, OverlaysPatchTheChangedTilesOnly) {
  const auto config = maps_.joined.config.get_child("mjolnir");
  const std::string dir = VALHALLA_BUILD_DIR "test/data/gurka_packages/split/overlays";
  fs::remove_all(dir);
  const auto report = PackageSet::WriteOverlays(config, dir);
  ASSERT_EQ(report.overlays.size(), 2u);
  EXPECT_EQ(report.stats.lost, 0u);

  // the runtime join, which rewrites every tile in memory
  const auto runtime = PackageSet::FromConfig(config);
  ASSERT_FALSE(runtime->uses_overlays());
  const auto joined = LoadAll(*runtime);

  size_t overlay_tiles = 0;
  for (const auto& package : maps_.packages) {
    // the tiles the join changed, compared with the package's own
    size_t changed = 0;
    for (const auto& [real, bytes] : TarTileBytes(package.tar)) {
      const auto copy = CopyOf(*runtime, package.name, real);
      ASSERT_TRUE(copy.is_valid()) << package.name << " " << real;
      changed += joined.at(copy) != bytes;
    }
    EXPECT_GT(changed, 0u) << package.name;
    overlay_tiles += changed;
    const auto written = std::find_if(report.overlays.begin(), report.overlays.end(),
                                      [&](const auto& o) { return o.name == package.name; });
    ASSERT_NE(written, report.overlays.end());
    EXPECT_EQ(written->tiles, changed) << package.name;
    EXPECT_EQ(written->bytes, fs::file_size(OverlayFile(dir, package.name))) << package.name;
    EXPECT_LT(written->bytes, fs::file_size(package.tar)) << package.name;
  }

  // with the overlays, loading a tile rewrites nothing and gives the same tiles
  const auto with_overlays = WithOverlays(maps_.joined, dir);
  const auto served = PackageSet::FromConfig(with_overlays.config.get_child("mjolnir"));
  ASSERT_TRUE(served->uses_overlays());
  EXPECT_EQ(served->join_key(), report.key);
  EXPECT_TRUE(LoadAll(*served) == joined);
  const auto stats = served->stats();
  EXPECT_EQ(stats.full_rewrites, 0u);
  EXPECT_EQ(stats.id_rewrites, 0u);
  EXPECT_EQ(stats.overlay_tiles, overlay_tiles);
  ExpectSingleGraphRoutes(with_overlays, maps_.single);
  GraphReader reader(with_overlays.config.get_child("mjolnir"));
  ExpectConsistentGraph(reader);
}

// Joined tiles are byte-identical on every platform. The grid of SplitJoin, shrunk into one
// level-2 tile: tile builds number the nodes of levels 0 and 1 in the order they visit the level-2
// tiles, which varies from build to build when there are several. The overlay files themselves
// aren't hashed: their patches are deflated, and zlib builds differ in what they write.
TEST(GoldenJoin, OverlayHashes) {
  const auto layout = gurka::detail::map_to_coordinates(SplitAscii(), 50, {5.2, 52.1});
  const double line = (layout.at("D").lng() + layout.at("E").lng()) / 2;
  const auto maps = BuildSplit(VALHALLA_BUILD_DIR "test/data/gurka_packages/golden", layout,
                               SplitWays(), {}, line, 500);
  const std::string dir = VALHALLA_BUILD_DIR "test/data/gurka_packages/golden/overlays";
  const auto report = PackageSet::WriteOverlays(maps.joined.config.get_child("mjolnir"), dir);
  EXPECT_EQ(report.stats.lost, 0u);
  // When the join or the tile format changes on purpose, update the hashes from this output.
  EXPECT_EQ(Fnv1a(ReadFile(maps.packages[0].tar)), 535886572131682744ull) << "west tiles";
  EXPECT_EQ(Fnv1a(ReadFile(maps.packages[1].tar)), 9302745409726607298ull) << "east tiles";
  // the tiles each package serves from its overlay, in tile order
  const auto served =
      PackageSet::FromConfig(WithOverlays(maps.joined, dir).config.get_child("mjolnir"));
  ASSERT_TRUE(served->uses_overlays());
  const auto tiles = LoadAll(*served);
  const auto joined_hash = [&](const Package& package) {
    std::string bytes;
    for (const auto& [real, original] : TarTileBytes(package.tar)) {
      bytes += tiles.at(CopyOf(*served, package.name, real));
    }
    return Fnv1a(bytes);
  };
  const auto patches_hash = [&](const std::string& name) {
    const auto overlay = std::find_if(report.overlays.begin(), report.overlays.end(),
                                      [&](const auto& o) { return o.name == name; });
    return overlay == report.overlays.end() ? 0 : overlay->content_hash;
  };
  EXPECT_EQ(patches_hash("west"), 10945871183313045016ull) << "west patches";
  EXPECT_EQ(patches_hash("east"), 3877405390329665803ull) << "east patches";
  EXPECT_EQ(joined_hash(maps.packages[0]), 3014461116394766175ull) << "west joined tiles";
  EXPECT_EQ(joined_hash(maps.packages[1]), 18440462409556805356ull) << "east joined tiles";
}

TEST_F(SplitJoin, MissingOrStaleOverlaysFallBackToTheRuntimeJoin) {
  const std::string dir = VALHALLA_BUILD_DIR "test/data/gurka_packages/split/stale";
  const auto config = maps_.joined.config.get_child("mjolnir");
  const auto runtime = LoadAll(*PackageSet::FromConfig(config));
  const auto expect_fallback = [&](const gurka::map& map, const std::string& what) {
    const auto set = PackageSet::FromConfig(map.config.get_child("mjolnir"));
    EXPECT_FALSE(set->uses_overlays()) << what;
    EXPECT_EQ(set->stats().overlay_tiles, 0u) << what;
    ExpectSingleGraphRoutes(map, maps_.single);
  };
  const auto rewrite = [&]() {
    fs::remove_all(dir);
    PackageSet::WriteOverlays(config, dir);
  };

  rewrite();
  ASSERT_TRUE(PackageSet::FromConfig(WithOverlays(maps_.joined, dir).config.get_child("mjolnir"))
                  ->uses_overlays());

  fs::remove(OverlayFile(dir, "east"));
  expect_fallback(WithOverlays(maps_.joined, dir), "missing overlay");
  EXPECT_TRUE(LoadAll(*PackageSet::FromConfig(
                  WithOverlays(maps_.joined, dir).config.get_child("mjolnir"))) == runtime);

  expect_fallback(WithOverlays(maps_.joined, dir + "/nonexistent"), "missing directory");

  // a changed package makes the overlays stale
  rewrite();
  auto newer = WithOverlays(maps_.joined, dir);
  newer.config.get_child("mjolnir.packages").back().second.put("build_time", 1);
  expect_fallback(newer, "changed build time");
  auto tolerance = WithOverlays(maps_.joined, dir);
  tolerance.config.put("mjolnir.package_join_tolerance", kJoinTolerance / 2);
  expect_fallback(tolerance, "changed join tolerance");

  // overlays of two different joins
  const auto other = dir + "_other";
  fs::remove_all(other);
  PackageSet::WriteOverlays(tolerance.config.get_child("mjolnir"), other);
  fs::copy_file(OverlayFile(other, "east"), OverlayFile(dir, "east"),
                fs::copy_options::overwrite_existing);
  expect_fallback(WithOverlays(maps_.joined, dir), "overlay of another join");

  // damaged overlays
  rewrite();
  const auto east = ReadFile(OverlayFile(dir, "east"));
  WriteFile(OverlayFile(dir, "east"), east.substr(0, east.size() / 2));
  expect_fallback(WithOverlays(maps_.joined, dir), "truncated overlay");
  WriteFile(OverlayFile(dir, "east"), std::string(east.size(), 'x'));
  expect_fallback(WithOverlays(maps_.joined, dir), "garbage overlay");
  WriteFile(OverlayFile(dir, "east"), east.substr(0, east.size() - 1));
  expect_fallback(WithOverlays(maps_.joined, dir), "overlay without its last byte");
  // one flipped bit in the middle of the patches
  auto flipped = east;
  flipped[flipped.size() * 3 / 4] ^= 1;
  WriteFile(OverlayFile(dir, "east"), flipped);
  expect_fallback(WithOverlays(maps_.joined, dir), "overlay with a flipped bit");
  WriteFile(OverlayFile(dir, "east"), east + std::string(1, '\0'));
  expect_fallback(WithOverlays(maps_.joined, dir), "overlay with a trailing byte");
}

TEST_F(SplitJoin, OverlaysDependOnlyOnThePackageSet) {
  const auto base = std::string(VALHALLA_BUILD_DIR "test/data/gurka_packages/split/identity");
  fs::remove_all(base);
  const auto config = maps_.joined.config.get_child("mjolnir");
  const auto first = PackageSet::WriteOverlays(config, base + "/first");

  // the packages in the other order
  auto reversed = config;
  boost::property_tree::ptree list;
  for (auto it = config.get_child("packages").rbegin(); it != config.get_child("packages").rend();
       ++it) {
    list.push_back(*it);
  }
  reversed.put_child("packages", list);
  const auto second = PackageSet::WriteOverlays(reversed, base + "/reversed");
  EXPECT_EQ(second.key, first.key);

  // another package set loaded earlier in the process, and still alive
  auto renamed = config;
  for (auto& entry : renamed.get_child("packages")) {
    entry.second.put("name", "other-" + entry.second.get<std::string>("name"));
    entry.second.put("build_time", 7);
  }
  GraphReader other_reader(renamed);
  const auto other = PackageSet::FromConfig(renamed);
  LoadAll(*other);
  const auto third = PackageSet::WriteOverlays(config, base + "/later");
  EXPECT_EQ(third.key, first.key);

  for (const auto& package : maps_.packages) {
    const auto expected = ReadFile(OverlayFile(base + "/first", package.name));
    EXPECT_FALSE(expected.empty());
    EXPECT_TRUE(ReadFile(OverlayFile(base + "/reversed", package.name)) == expected) << package.name;
    EXPECT_TRUE(ReadFile(OverlayFile(base + "/later", package.name)) == expected) << package.name;
  }
  // and the overlays serve a set loaded next to the other one
  const auto served =
      PackageSet::FromConfig(WithOverlays(maps_.joined, base + "/later").config.get_child("mjolnir"));
  EXPECT_TRUE(served->uses_overlays());
}

// The overlays are written by several threads, in tile order: the thread count changes how long a
// join takes and nothing it writes.
TEST_F(SplitJoin, OverlaysDoNotDependOnTheThreadCount) {
  const auto base = std::string(VALHALLA_BUILD_DIR "test/data/gurka_packages/split/threads");
  fs::remove_all(base);
  auto config = maps_.joined.config.get_child("mjolnir");
  std::map<size_t, PackageSet::JoinReport> reports;
  for (const size_t threads : {1u, 2u, 7u, 64u}) {
    config.put("concurrency", threads);
    reports[threads] = PackageSet::WriteOverlays(config, base + "/" + std::to_string(threads));
  }
  for (const auto& [threads, report] : reports) {
    EXPECT_EQ(report.key, reports.at(1).key) << threads << " threads";
    EXPECT_EQ(report.stats.joined, reports.at(1).stats.joined) << threads << " threads";
    EXPECT_EQ(report.stats.lost, reports.at(1).stats.lost) << threads << " threads";
    for (const auto& package : maps_.packages) {
      const auto expected = ReadFile(OverlayFile(base + "/1", package.name));
      EXPECT_FALSE(expected.empty());
      EXPECT_TRUE(ReadFile(OverlayFile(base + "/" + std::to_string(threads), package.name)) ==
                  expected)
          << package.name << " with " << threads << " threads";
    }
  }
  // the golden overlays are the ones a multi-threaded join writes too
  const auto served =
      PackageSet::FromConfig(WithOverlays(maps_.joined, base + "/7").config.get_child("mjolnir"));
  EXPECT_TRUE(served->uses_overlays());
}

TEST_F(SplitJoin, TileRefsMustMatchTheirTar) {
  const auto& west = maps_.packages[0];
  const auto& east = maps_.packages[1];
  const auto refs = ReadFile(east.refs);
  ASSERT_GE(refs.size(), 16u);
  // magic and format version
  EXPECT_EQ(refs.substr(0, 4), "VTRF");
  uint32_t version = 0;
  std::memcpy(&version, refs.data() + 4, sizeof(version));
  EXPECT_EQ(version, 2u);

  const auto config = maps_.joined.config.get_child("mjolnir");
  const auto set = PackageSet::FromConfig(config);
  const auto expected = LoadAll(*set);
  const auto baseline = set->stats().full_rewrites;
  ASSERT_LT(baseline, expected.size()) << "the index saves work";
  size_t east_tiles = 0;
  for (const auto& id : set->AllTiles()) {
    east_tiles += set->name(set->PackageOf(id)) == "east";
  }

  const std::string dir = VALHALLA_BUILD_DIR "test/data/gurka_packages/split/";
  const auto with_refs = [&](const std::string& file) {
    auto pt = config;
    for (auto& entry : pt.get_child("packages")) {
      if (entry.second.get<std::string>("name") == "east") {
        entry.second.put("tile_refs", file);
      }
    }
    return pt;
  };
  // every east tile gets a full rewrite, and the tiles come out the same
  const auto ignored = [&](const boost::property_tree::ptree& pt, const std::string& what) {
    const auto set = PackageSet::FromConfig(pt);
    EXPECT_TRUE(LoadAll(*set) == expected) << what;
    EXPECT_GE(set->stats().full_rewrites, east_tiles) << what;
    EXPECT_GT(set->stats().full_rewrites, baseline) << what;
  };
  ignored(with_refs(west.refs), "another package's index");
  WriteFile(dir + "old_format.bin", refs.substr(16));
  ignored(with_refs(dir + "old_format.bin"), "index without header");
  WriteFile(dir + "garbage.bin", std::string(refs.size(), 'x'));
  ignored(with_refs(dir + "garbage.bin"), "garbage");
  WriteFile(dir + "truncated.bin", refs.substr(0, refs.size() - 2));
  ignored(with_refs(dir + "truncated.bin"), "truncated index");
  ignored(with_refs(dir + "nonexistent.bin"), "missing index");
}

// ---------------------------------------------------------------------------------------------
// Two crossings; east's copy of C lies 15 m north of west's, more than the join tolerance.

TEST(MovedTwin, CrossingIsDisabledNotLinked) {
  const std::string ascii = R"(
    A-----B---C-----D
    |     |   |     |
    E-----F---G-----H
  )";
  const auto layout = gurka::detail::map_to_coordinates(ascii, 100, {5.1, 52.1});
  auto moved = layout;
  moved["C"] = midgard::PointLL(layout.at("C").lng(), layout.at("C").lat() + 15 / kMetersPerDegLat);
  const gurka::ways ways = {
      {"AB", {{"highway", "residential"}}}, {"BC", {{"highway", "residential"}}},
      {"CD", {{"highway", "residential"}}}, {"EF", {{"highway", "residential"}}},
      {"FG", {{"highway", "residential"}}}, {"GH", {{"highway", "residential"}}},
      {"AE", {{"highway", "residential"}}}, {"BF", {{"highway", "residential"}}},
      {"CG", {{"highway", "residential"}}}, {"DH", {{"highway", "residential"}}},
  };
  const double line = (layout.at("B").lng() + layout.at("C").lng()) / 2;
  const auto maps = BuildSplit(VALHALLA_BUILD_DIR "test/data/gurka_packages/moved", layout, ways, {},
                               line, 3000, moved);

  const auto set = PackageSet::FromConfig(maps.joined.config.get_child("mjolnir"));
  for (const auto& id : set->AllTiles()) {
    set->LoadTile(id);
  }
  EXPECT_EQ(set->stats().lost, 2u) << "B-C in both directions";
  EXPECT_EQ(set->stats().joined, 2u) << "F-G in both directions";

  for (const auto& [from, to] : {std::make_pair("A", "D"), std::make_pair("D", "A")}) {
    const auto route = RouteWays(maps.joined, from, to);
    EXPECT_EQ(std::count(route.begin(), route.end(), WayId(maps.ways, "BC")), 0)
        << from << " -> " << to;
    EXPECT_EQ(std::count(route.begin(), route.end(), WayId(maps.ways, "FG")), 1)
        << from << " -> " << to;
  }
  GraphReader reader(maps.joined.config.get_child("mjolnir"));
  ExpectConsistentGraph(reader);
}

// ---------------------------------------------------------------------------------------------
// Complex restrictions (via a way) across the join. The via way X-Y crosses the line, so the from
// and to edges lie in different packages; one restriction per direction, so the from edge is in
// west for one and in east for the other. The south road C-Z-W-D is the detour.

struct RestrictionCase {
  std::string name;
  // the maneuver each restriction forbids, and the restriction tags
  std::string restriction_east, restriction_west;
};

class CrossingRestriction : public ::testing::TestWithParam<RestrictionCase> {};

TEST_P(CrossingRestriction, IsEnforced) {
  const std::string ascii = R"(
    M-N                    A----X-----Y----B                    O-P
                                |     |
                           C----Z-----W----D
  )";
  const auto layout = gurka::detail::map_to_coordinates(ascii, 100, {5.1, 52.1});
  const gurka::ways ways = {
      {"MN", {{"highway", "residential"}}}, {"OP", {{"highway", "residential"}}},
      {"AX", {{"highway", "residential"}}}, {"XY", {{"highway", "residential"}}},
      {"YB", {{"highway", "residential"}}}, {"CZ", {{"highway", "residential"}}},
      {"ZW", {{"highway", "residential"}}}, {"WD", {{"highway", "residential"}}},
      {"XZ", {{"highway", "residential"}}}, {"YW", {{"highway", "residential"}}},
  };
  const auto& param = GetParam();
  // eastbound: from A-X via X-Y; westbound: from B-Y via Y-X
  const bool only = param.name == "only";
  const gurka::relations relations = {
      {{{gurka::way_member, "AX", "from"},
        {gurka::way_member, "XY", "via"},
        {gurka::way_member, only ? "YW" : "YB", "to"}},
       {{"type", "restriction"}, {"restriction", param.restriction_east}}},
      {{{gurka::way_member, "YB", "from"},
        {gurka::way_member, "XY", "via"},
        {gurka::way_member, only ? "XZ" : "AX", "to"}},
       {{"type", "restriction"}, {"restriction", param.restriction_west}}},
  };
  const double line = (layout.at("X").lng() + layout.at("Y").lng()) / 2;
  const auto maps =
      BuildSplit(VALHALLA_BUILD_DIR "test/data/gurka_packages/restriction_" + param.name, layout,
                 ways, relations, line, 1500);

  const auto ax = WayId(maps.ways, "AX"), xy = WayId(maps.ways, "XY"), yb = WayId(maps.ways, "YB");
  const auto contains = [](const std::vector<uint64_t>& route, std::vector<uint64_t> maneuver) {
    return std::search(route.begin(), route.end(), maneuver.begin(), maneuver.end()) != route.end();
  };
  // the single graph enforces them
  EXPECT_FALSE(contains(RouteWays(maps.single, "A", "B"), {ax, xy, yb}));
  EXPECT_FALSE(contains(RouteWays(maps.single, "B", "A"), {yb, xy, ax}));
  // and so do the packages, with the from edge in west, then in east
  const auto east = RouteWays(maps.joined, "A", "B");
  const auto west = RouteWays(maps.joined, "B", "A");
  EXPECT_FALSE(contains(east, {ax, xy, yb}));
  EXPECT_FALSE(contains(west, {yb, xy, ax}));
  EXPECT_EQ(east, RouteWays(maps.single, "A", "B"));
  EXPECT_EQ(west, RouteWays(maps.single, "B", "A"));
  // without the maneuver in question the via way is still usable across the join
  EXPECT_TRUE(contains(RouteWays(maps.joined, "A", "W"), {ax, xy}));
}

INSTANTIATE_TEST_SUITE_P(Restrictions,
                         CrossingRestriction,
                         ::testing::Values(RestrictionCase{"no", "no_straight_on", "no_straight_on"},
                                           RestrictionCase{"only", "only_right_turn",
                                                           "only_left_turn"}),
                         [](const auto& info) { return info.param.name; });

// ---------------------------------------------------------------------------------------------
// Six packages along one road, all in one level-2 tile.

TEST(SixPackages, HoldingOneTileLoadAndRoute) {
  const std::string ascii = R"(
    A--B--C--D--E--F--G--H--I--J--K--L--M
  )";
  const auto layout = gurka::detail::map_to_coordinates(ascii, 100, {5.1, 52.1});
  const std::string nodes = "ABCDEFGHIJKLM";
  gurka::ways ways;
  for (size_t i = 0; i + 1 < nodes.size(); ++i) {
    ways[nodes.substr(i, 2)] = {{"highway", "residential"}};
  }
  ways = WithWayIds(ways);
  const std::string dir = VALHALLA_BUILD_DIR "test/data/gurka_packages/six";
  fs::remove_all(dir);
  fs::create_directories(dir);
  // package k owns the nodes 2k and 2k + 1 and is built from the ways within one node of them
  std::vector<Package> packages;
  const auto lng = [&](size_t i) { return layout.at(std::string(1, nodes[i])).lng(); };
  const double spacing = lng(1) - lng(0);
  for (size_t k = 0; k < 6; ++k) {
    const double min = lng(2 * k) - spacing / 2, max = lng(2 * k + 1) + spacing / 2;
    packages.push_back(BuildPackage(dir, "p" + std::to_string(k), layout,
                                    WaysWithin(ways, layout, min - spacing, max + spacing), {},
                                    {k == 0 ? 4.0 : min, 51.0, k == 5 ? 6.0 : max, 53.0}));
  }
  const auto map = PackageMap(packages, dir, layout);

  GraphReader reader(map.config.get_child("mjolnir"));
  const auto tile = TileHierarchy::GetGraphId(layout.at("A"), 2);
  EXPECT_EQ(reader.GetTileCopies(tile).size(), 6u);
  const auto set = PackageSet::FromConfig(map.config.get_child("mjolnir"));
  for (const auto& id : set->AllTiles()) {
    set->LoadTile(id);
  }
  EXPECT_EQ(set->stats().lost, 0u);
  EXPECT_EQ(set->stats().joined, 10u) << "five package borders, both directions";
  EXPECT_NEAR(RouteKm(map, "A", "M"), Km(layout, "A", "M"), 0.01);
  EXPECT_NEAR(RouteKm(map, "M", "A"), Km(layout, "M", "A"), 0.01);
  ExpectConsistentGraph(reader);
}

// ---------------------------------------------------------------------------------------------
// Three packages, west, middle, and east, along one road, all in one level-2 tile. Each owns a run
// of nodes and is built from the ways within one node of its box, so west and east reach into
// middle's box without meeting each other.

class ThreePackages : public ::testing::Test {
protected:
  static void SetUpTestSuite() {
    const std::string ascii = R"(
      A--B--C--D--E--F--G--H--I--J--K--L--M
    )";
    layout_ = gurka::detail::map_to_coordinates(ascii, 100, {5.1, 52.1});
    for (size_t i = 0; i + 1 < kNodes.size(); ++i) {
      ways_[kNodes.substr(i, 2)] = {{"highway", "residential"}};
    }
    ways_ = WithWayIds(ways_);
    const std::string dir = VALHALLA_BUILD_DIR "test/data/gurka_packages/three";
    fs::remove_all(dir);
    fs::create_directories(dir);
    single_ = gurka::buildtiles(layout_, ways_, {}, {}, dir + "/single");
    const double spacing = Lng(1) - Lng(0);
    const auto build = [&](const std::string& name, size_t first, size_t last, double min_lng,
                           double max_lng) {
      const double min = Lng(first) - spacing / 2, max = Lng(last) + spacing / 2;
      return BuildPackage(dir, name, layout_,
                          WaysWithin(ways_, layout_, min - spacing, max + spacing), {},
                          {min_lng ? min_lng : min, 51.0, max_lng ? max_lng : max, 53.0});
    };
    const auto west = build("west", 0, 3, 4.0, 0);
    const auto middle = build("middle", 4, 8, 0, 0);
    const auto east = build("east", 9, 12, 0, 6.0);
    without_middle_ = PackageMap({west, east}, dir + "/without_middle", layout_);
    all_ = PackageMap({west, middle, east}, dir + "/all", layout_);
  }

  static double Lng(size_t i) {
    return layout_.at(std::string(1, kNodes[i])).lng();
  }

  static const std::string kNodes;
  static gurka::nodelayout layout_;
  static gurka::ways ways_;
  static gurka::map single_, without_middle_, all_;
};
const std::string ThreePackages::kNodes = "ABCDEFGHIJKLM";
gurka::nodelayout ThreePackages::layout_;
gurka::ways ThreePackages::ways_;
gurka::map ThreePackages::single_, ThreePackages::without_middle_, ThreePackages::all_;

TEST_F(ThreePackages, DeletingOneEqualsASetThatNeverHadIt) {
  // the same routes with all three packages installed, so the test can tell the sets apart
  for (const auto& map : {without_middle_, all_}) {
    for (const auto& [from, to] : {std::make_pair("A", "D"), std::make_pair("D", "A"),
                                   std::make_pair("J", "M"), std::make_pair("M", "J")}) {
      const auto expected = Route(single_, from, to);
      const auto actual = Route(map, from, to);
      EXPECT_EQ(actual.edges, expected.edges) << from << " -> " << to;
      EXPECT_NEAR(actual.seconds, expected.seconds, 1e-3) << from << " -> " << to;
    }
  }
  GraphReader reader(without_middle_.config.get_child("mjolnir"));
  ExpectConsistentGraph(reader);
}

TEST_F(ThreePackages, BorderRoadsEndWithTheRemainingData) {
  // F is west's last node and lies in middle's box
  EXPECT_NEAR(RouteKm(without_middle_, "A", "F"), Km(layout_, "A", "F"), 0.01);
  EXPECT_NEAR(RouteKm(without_middle_, "F", "A"), Km(layout_, "F", "A"), 0.01);
  EXPECT_THROW(RouteKm(without_middle_, "A", "M"), std::exception);
  EXPECT_THROW(RouteKm(without_middle_, "M", "A"), std::exception);
  // with middle installed, the same trip goes through
  EXPECT_NEAR(RouteKm(all_, "A", "M"), Km(layout_, "A", "M"), 0.01);
}

// ---------------------------------------------------------------------------------------------
// An outer package around the whole map and a newer inner one inside it, like a country and a
// region of it.

TEST(NestedPackages, DeletingTheInnerGivesTheAreaBack) {
  const std::string ascii = R"(
    A--B--C--D--E--F--G
  )";
  const auto layout = gurka::detail::map_to_coordinates(ascii, 100, {5.1, 52.1});
  gurka::ways ways;
  const std::string nodes = "ABCDEFG";
  for (size_t i = 0; i + 1 < nodes.size(); ++i) {
    ways[nodes.substr(i, 2)] = {{"highway", "residential"}};
  }
  ways = WithWayIds(ways);
  const std::string dir = VALHALLA_BUILD_DIR "test/data/gurka_packages/nested";
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto single = gurka::buildtiles(layout, ways, {}, {}, dir + "/single");

  const double spacing = layout.at("B").lng() - layout.at("A").lng();
  const double min = layout.at("C").lng() - spacing / 2, max = layout.at("E").lng() + spacing / 2;
  const auto outer = BuildPackage(dir, "outer", layout, ways, {}, {4.0, 51.0, 6.0, 53.0}, 1);
  const auto inner =
      BuildPackage(dir, "inner", layout, WaysWithin(ways, layout, min - spacing, max + spacing), {},
                   {min, 51.0, max, 53.0}, 2);

  const auto both = PackageMap({outer, inner}, dir + "/both", layout);
  const auto set = PackageSet::FromConfig(both.config.get_child("mjolnir"));
  const auto owner = [&](const std::string& node) { return set->name(set->Owner(layout.at(node))); };
  EXPECT_EQ(owner("D"), "inner");
  EXPECT_EQ(owner("A"), "outer");

  const auto outer_alone = PackageMap({outer}, dir + "/outer_alone", layout);
  const auto expected = Route(single, "C", "E");
  const auto actual = Route(outer_alone, "C", "E");
  EXPECT_EQ(actual.edges, expected.edges);
  EXPECT_NEAR(actual.seconds, expected.seconds, 1e-3);
}

// ---------------------------------------------------------------------------------------------
// Twenty packages whose polygons all contain one point.

class OverlappingPackages : public ::testing::Test {
protected:
  static void SetUpTestSuite() {
    dir_ = VALHALLA_BUILD_DIR "test/data/gurka_packages/overlap";
    fs::remove_all(dir_);
    fs::create_directories(dir_);
    const auto layout = gurka::detail::map_to_coordinates("A---B", 100, kCenter);
    base_ = BuildPackage(dir_, "base", layout, {{"AB", {{"highway", "residential"}}}}, {},
                         {5.0, 52.0, 5.1, 52.1});
  }

  // Twenty packages over the same tiles, named in index order; package i has a square polygon
  // around the center that is i + 1 hundredths of a degree wide, and the given build time.
  static boost::property_tree::ptree Config(const std::function<int64_t(size_t)>& build_time) {
    std::vector<Package> packages;
    for (size_t i = 0; i < 20; ++i) {
      auto package = base_;
      package.name = (i < 10 ? "p0" : "p") + std::to_string(i);
      package.poly = dir_ + "/" + package.name + ".poly";
      const double half = 0.005 * (i + 1);
      WriteBox(package.poly, kCenter.lng() - half, kCenter.lat() - half, kCenter.lng() + half,
               kCenter.lat() + half);
      package.build_time = build_time(i);
      packages.push_back(package);
    }
    boost::property_tree::ptree config;
    config.put_child("packages", PackageList(packages));
    config.put("package_join_tolerance", kJoinTolerance);
    return config;
  }

  static const midgard::PointLL kCenter;
  static std::string dir_;
  static Package base_;
};
const midgard::PointLL OverlappingPackages::kCenter{5.33, 52.33};
std::string OverlappingPackages::dir_;
Package OverlappingPackages::base_;

TEST_F(OverlappingPackages, DeepestWinsAmongEqualBuildTimes) {
  const auto set = PackageSet::FromConfig(Config([](size_t) { return 100; }));
  // the widest polygon contains the center most deeply
  EXPECT_EQ(set->Owner(kCenter), 19);
  EXPECT_EQ(set->Owner({kCenter.lng() + 0.001, kCenter.lat() - 0.002}), 19);
  EXPECT_EQ(set->Copies(TileHierarchy::GetGraphId(kCenter, 2)).size(), 20u);
}

TEST_F(OverlappingPackages, NewestWins) {
  const auto set = PackageSet::FromConfig(Config([](size_t i) { return i == 17 ? 200 : 100; }));
  EXPECT_EQ(set->Owner(kCenter), 17);
  // outside its polygon the deepest of the others owns
  EXPECT_EQ(set->Owner({kCenter.lng() + 0.093, kCenter.lat()}), 19);
}

// ---------------------------------------------------------------------------------------------
// A package whose tar lacks tiles its own tiles point into. Package "cut" ships only the level-2
// tile of A and B; its build also had the level-2 tile of C and D (east of a tile boundary) and
// the level-0 tile of the primary road. Package "full" ships every tile, so the ids "cut" points
// with are ids of real tiles in the set, just not of "cut"'s copies.

TEST(MissingTiles, ReferencesIntoThemAreDisabled) {
  const std::string ascii = R"(
    A----B-------C----D
         |
         E
  )";
  // a level-2 tile boundary runs between B and C
  const double boundary = 5.25;
  const auto probe = gurka::detail::map_to_coordinates(ascii, 100, {5.0, 52.1});
  const double shift = 5.0 + boundary - (probe.at("B").lng() + probe.at("C").lng()) / 2;
  const auto layout = gurka::detail::map_to_coordinates(ascii, 100, {shift, 52.1});
  // B is on level 2 and, for the primary road, on level 0: it has transitions between them
  const auto ways = WithWayIds({
      {"AB", {{"highway", "residential"}}},
      {"BC", {{"highway", "residential"}}},
      {"CD", {{"highway", "residential"}}},
      {"BE", {{"highway", "primary"}}},
  });
  const gurka::relations relations = {
      {{{gurka::way_member, "CD", "from"},
        {gurka::way_member, "BC", "via"},
        {gurka::way_member, "AB", "to"}},
       {{"type", "restriction"}, {"restriction", "no_straight_on"}}},
  };
  const auto tile_a = TileHierarchy::GetGraphId(layout.at("A"), 2);
  const auto tile_c = TileHierarchy::GetGraphId(layout.at("C"), 2);
  ASSERT_NE(tile_a, tile_c);
  const auto level0 = TileHierarchy::GetGraphId(layout.at("A"), 0);

  const std::string dir = VALHALLA_BUILD_DIR "test/data/gurka_packages/missing";
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto full = BuildPackage(dir, "full", layout, ways, relations, {4.0, 51.0, 6.0, 53.0}, 1);
  const auto cut =
      BuildPackage(dir, "cut", layout, ways, relations, {4.9, 51.9, boundary - 0.001, 52.3}, 2,
                   {tile_c, level0, TileHierarchy::GetGraphId(layout.at("A"), 1)});
  const auto map = PackageMap({full, cut}, dir, layout);

  GraphReader reader(map.config.get_child("mjolnir"));
  const auto& set = *reader.package_set();
  const auto copy = CopyOf(set, "cut", tile_a);
  ASSERT_TRUE(copy.is_valid());
  const auto tile = reader.GetGraphTile(copy);
  ASSERT_TRUE(tile);
  const int cut_index = set.PackageOf(copy);

  // a reference is either into one of cut's copies, or into no tile at all
  const auto expect_own_or_none = [&](const GraphId& id, const std::string& what) {
    if (id.is_valid() && reader.DoesTileExist(id)) {
      EXPECT_EQ(set.PackageOf(id.tile_base()), cut_index) << what << " " << id;
    }
  };
  size_t transitions = 0, restrictions = 0;
  for (uint32_t t = 0; t < tile->header()->transitioncount(); ++t) {
    expect_own_or_none(tile->transition(t)->endnode(), "transition");
    ++transitions;
  }
  EXPECT_GT(transitions, 0u);
  for (size_t b = 0; b < kBinCount; ++b) {
    for (const auto& id : tile->GetBin(b)) {
      expect_own_or_none(id, "bin entry");
    }
  }
  for (uint32_t e = 0; e < tile->header()->directededgecount(); ++e) {
    const auto* edge = tile->directededge(e);
    if (edge->forwardaccess() || edge->reverseaccess()) {
      EXPECT_TRUE(reader.DoesTileExist(edge->endnode())) << "usable edge " << e;
    }
    for (const bool forward : {true, false}) {
      for (const auto& restriction :
           tile->GetComplexRestrictions(forward, GraphId(copy.tileid(), copy.level(), e),
                                        kAllAccess)) {
        expect_own_or_none(restriction.from_graphid(), "restriction from");
        expect_own_or_none(restriction.to_graphid(), "restriction to");
        restriction.WalkVias([&](const GraphId* via) {
          expect_own_or_none(*via, "restriction via");
          return WalkingVia::KeepWalking;
        });
        ++restrictions;
      }
    }
  }
  EXPECT_GT(restrictions, 0u);
  // routing inside cut's tile still works
  EXPECT_NEAR(RouteKm(map, "A", "B"), Km(layout, "A", "B"), 0.01);
}

} // namespace
