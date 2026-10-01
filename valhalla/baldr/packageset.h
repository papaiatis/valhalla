#pragma once

#include <valhalla/baldr/graphid.h>
#include <valhalla/baldr/graphtileptr.h>
#include <valhalla/midgard/aabb2.h>
#include <valhalla/midgard/pointll.h>

#include <boost/property_tree/ptree_fwd.hpp>

#include <atomic>
#include <cstdint>
#include <limits>
#include <memory>
#include <shared_mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace valhalla {
namespace baldr {

class GraphTile;
class DirectedEdge;

/**
 * Region polygon read from an osmosis .poly file. Answers point containment and distance to the
 * boundary. A file may hold several rings; holes (rings whose name starts with '!') are handled
 * by even-odd parity. Throws std::runtime_error on a malformed file.
 */
class RegionPolygon {
public:
  explicit RegionPolygon(const std::string& poly_file);

  bool Contains(const midgard::PointLL& p) const;

  // Distance in meters from p to the nearest boundary segment.
  double Distance(const midgard::PointLL& p) const;

  const midgard::AABB2<midgard::PointLL>& bbox() const {
    return bbox_;
  }

private:
  struct Segment {
    midgard::PointLL a, b;
  };
  std::vector<Segment> segments_;
  midgard::AABB2<midgard::PointLL> bbox_;

  // Segment index: square cells of kCellDeg degrees -> segment indices.
  std::unordered_map<int64_t, std::vector<uint32_t>> cells_;
  // Row bands of kCellDeg degrees -> segment indices overlapping the band (for ray casting).
  std::unordered_map<int32_t, std::vector<uint32_t>> rows_;
  // Range of the cells that hold segments.
  int32_t row_min_ = std::numeric_limits<int32_t>::max();
  int32_t row_max_ = std::numeric_limits<int32_t>::min();
  int32_t col_min_ = std::numeric_limits<int32_t>::max();
  int32_t col_max_ = std::numeric_limits<int32_t>::min();

  static int32_t Row(double lat);
  static int32_t Col(double lng);
  static int64_t Key(int32_t row, int32_t col) {
    return (static_cast<int64_t>(row) << 32) | static_cast<uint32_t>(col);
  }
};

/**
 * A set of routing packages, each built as an independent Valhalla graph from its own extract.
 *
 * Tiles of different packages occupy the same geographic tile ids but are numbered by different
 * builds. Every package copy of a tile gets its own tile id: the first copy ("slot 0") keeps the
 * real id, every further copy gets an unused id above the level's tile grid. A tile is loaded
 * from its package and all GraphIds it stores are rewritten to the ids of the same package's
 * copies, so every in-memory reference stays inside one build. A reference into a tile the
 * package does not ship points into no tile at all.
 *
 * Packages are joined through ownership: every location is owned by the newest package whose
 * polygon contains it; among packages with equal build times, by the one containing it most
 * deeply; outside every polygon, by the nearest package. A package's edges that start at a node it
 * does not own are unreachable: they are dropped from the edge bins, and the ones leading back to
 * an owned node are disabled. An edge that starts at an owned node and ends at a node owned by
 * another package is redirected to that package's twin node: the mutually nearest node within the
 * join tolerance whose edge back along the same way ends at the twin of our start node. Otherwise
 * the edge is disabled. Complex restrictions refer to the copies of their edges that routes use.
 *
 * The join happens when a tile is first loaded (the runtime join, rewriting the tile in a private
 * copy-on-write mapping of its package), or once ahead of time: WriteOverlays writes, per package,
 * the tiles the join changed into an overlay, and a set configured with those overlays serves
 * every tile as stored. Overlays that do not match the configured packages are ignored with a
 * warning, and the set falls back to the runtime join.
 *
 * Configuration, in the mjolnir section:
 *   packages                list of {name, tile_extract, polygon, build_time, tile_refs}:
 *                           a unique name, which also names the package's overlay; its tile tar;
 *                           its .poly file; its build time (seconds since the epoch), newer
 *                           builds own more; and optionally its tile reference index
 *                           (valhalla_build_tile_refs), which saves rewriting tiles
 *   package_join_tolerance  meters a node and its twin may lie apart (OSM edits move nodes)
 *   package_joined          directory holding the overlays <name>.joined (valhalla_join_packages)
 * Package files must be replaced by renaming, never rewritten in place: they stay mapped.
 *
 * One PackageSet serves every thread of the process; it needs a build with
 * ENABLE_THREAD_SAFE_TILE_REF_COUNT, as its tiles are shared between threads.
 */
class PackageSet {
public:
  /**
   * Builds the package set from the "packages" array of the mjolnir config and registers its copy
   * ids in the process. Returns nullptr when the config has no packages. Throws
   * std::runtime_error on missing or malformed package files.
   */
  static std::unique_ptr<PackageSet> FromConfig(const boost::property_tree::ptree& pt);

  struct Stats {
    uint64_t joined = 0;        // edges redirected to another package
    uint64_t lost = 0;          // crossing edges without an exact twin (disabled)
    uint64_t disabled = 0;      // edges starting at a node owned by another package
    uint64_t clean_tiles = 0;   // tiles served unmodified from their package
    uint64_t id_rewrites = 0;   // tiles with only their ids remapped
    uint64_t full_rewrites = 0; // tiles fully joined
    uint64_t overlay_tiles = 0; // tiles served from an overlay
  };

  struct Overlay {
    std::string name;   // package name
    uint64_t tiles = 0; // tiles the join changed
    uint64_t bytes = 0; // their size
  };
  struct JoinReport {
    uint64_t key = 0; // join key of the package set
    Stats stats;
    std::vector<Overlay> overlays;
  };

  /**
   * Joins every tile of the configured packages, one tile at a time, and writes each package's
   * overlay to dir/<name>.joined: a tar whose first entry holds the join key and the copy ids the
   * join used, followed by every tile the join changed. The overlays are written under temporary
   * names and renamed when all are complete. They depend only on the package set, not on the
   * order of the configuration or on other package sets in the process. Throws
   * std::runtime_error on failure, leaving the previous overlays in place.
   */
  static JoinReport WriteOverlays(const boost::property_tree::ptree& pt, const std::string& dir);

  ~PackageSet();
  PackageSet(const PackageSet&) = delete;
  PackageSet& operator=(const PackageSet&) = delete;

  // Loads the tile with the given (copy) tile base id, joining it if needed. nullptr if no such
  // copy. Throws std::runtime_error when the tile cannot be joined.
  graph_tile_ptr LoadTile(const GraphId& tile_base) const;

  // Tile base ids of every package copy of the given tile (real or copy id on input).
  std::vector<GraphId> Copies(const GraphId& tile) const;

  // Tile base ids of every copy of every tile of every package, in id order.
  std::vector<GraphId> AllTiles() const;

  bool Exists(const GraphId& tile_base) const;

  // Index of the package whose copy the tile is, or -1.
  int PackageOf(const GraphId& tile) const;

  // Packages are indexed in name order.
  size_t size() const;
  const std::string& name(size_t package) const;

  // Owning package index for a location (see the class comment).
  int Owner(const midgard::PointLL& p) const;

  // Whether tiles are served from overlays (see WriteOverlays).
  bool uses_overlays() const {
    return overlays_;
  }

  // Hash identifying the package set: every package's name, build time, polygon, and tiles, the
  // join tolerance, and the version of the join. Overlays with another key are stale.
  uint64_t join_key() const {
    return join_key_;
  }

  Stats stats() const;

  /**
   * Writes the tile reference index (the "tile_refs" file of a package) of a tile tar. It starts
   * with a magic number, its format version, and the fingerprint of the tar it indexes. Per tile
   * it holds the bounding box of the nodes outside the tile that the tile's GraphIds lead to (end
   * nodes of its edges, start nodes of the edges in its bins and complex restrictions) and every
   * other tile its GraphIds point into. With it, LoadTile serves tiles that need no rewrite
   * straight from the package. Returns the number of tile references written. Throws
   * std::runtime_error on failure.
   */
  static uint64_t WriteTileRefs(const std::string& tile_extract, const std::string& file);

  // The GraphId in the real tile of a GraphId in a package copy (unchanged for real ids).
  static GraphId Real(const GraphId& id);

  // Which copy of its real tile a GraphId is in: 0 for real ids.
  static uint32_t CopyIndex(const GraphId& id);

  static constexpr const char* kOverlaySuffix = ".joined";

private:
  struct Package;
  struct TileEntry;
  // The tile id of copy `copy` (> 0) of the real tile (level, tileid).
  struct CopyId {
    uint32_t level = 0, tileid = 0, copy = 0, id = 0;
    bool operator==(const CopyId& other) const {
      return level == other.level && tileid == other.tileid && copy == other.copy && id == other.id;
    }
  };
  PackageSet();

  // The process-wide table of copy ids, which GraphTile::BoundingBox and
  // TileHierarchy::GetGraphIdBoundingBox read through Real. Registering replaces the ids other live
  // sets hold for other tiles with free ones and returns false if it had to.
  static bool RegisterCopyIds(std::vector<CopyId>& ids);
  static void UnregisterCopyIds(const std::vector<CopyId>& ids);

  // The package set with the copy ids that depend only on the set; not registered in the process.
  static std::unique_ptr<PackageSet> Build(const boost::property_tree::ptree& pt);
  // Registers the copy ids in the process; ids another live set holds for other tiles are
  // replaced by free ones.
  void Register();
  // Indexes every copy by its tile id.
  void IndexCopies();
  // Serves tiles from the overlays in dir if they all match the set.
  void OpenOverlays(const std::string& dir);
  // Every copy id of the set beyond slot 0.
  std::vector<CopyId> CopyIds() const;

  // What loading a package's copy of a tile takes, from the bake-time tile reference index:
  //  kNone      the tile is self-clean (slot 0 and wholly owned by the package), the package owns
  //             every node its edges, bins, and restrictions lead to, and every tile it references
  //             is the package's own and in slot 0: serve it unmodified
  //  kLinksOnly as kNone, but a tile on another level is in a slot > 0: remap the header,
  //             transitions, bins, and complex restrictions; node and edge arrays are not read
  //             (directed edges only connect nodes of their own level)
  //  kIdsOnly   as kLinksOnly, and a tile on the same level is in a slot > 0: remap end nodes too
  //  kFull      anything else, or no index shipped
  enum class Work { kNone, kLinksOnly, kIdsOnly, kFull };

  // Id in the package's copy of `real` (a GraphId from package `pkg`'s own tiles); an id in no
  // tile if the package does not ship the tile.
  GraphId Virtual(uint32_t pkg, const GraphId& real) const;

  // Rewrites a tile of package pkg in place, doing the given work.
  void Rewrite(uint32_t pkg, GraphTile& tile, Work work) const;
  // Owners of the nodes of the tile being rewritten, found once per node.
  class NodeOwners;
  // The kFull part of Rewrite: joins or disables the edges of a tile with the given real id.
  void JoinEdges(uint32_t pkg, GraphTile& tile, const GraphId& real_base, NodeOwners& owners) const;
  // An edge id of a complex restriction of package pkg, as the copy routes use: the copy of the
  // package owning its start node.
  GraphId RestrictionEdge(uint32_t pkg, const GraphId& real) const;

  // A package's copy of a tile needs no rewrite by itself: slot 0 and wholly owned by the package.
  bool SelfClean(uint32_t pkg, const TileEntry& entry) const;
  Work ComputeWork(uint32_t pkg, const TileEntry& entry) const;
  // Whether no edge of a tile whose single owner is tile_owner can end at a node of the package,
  // from the reference index.
  bool ReachesOnlyForeign(uint32_t pkg, uint64_t tile_base, int tile_owner) const;

  struct Twin {
    GraphId endnode;    // node id in the other package's copy
    uint32_t opp_index; // index of the opposing edge at that node
    uint32_t opp_local_idx;
    uint32_t restrictions; // simple restriction mask of the other package's copy of the edge
    uint32_t reverseaccess;
    const char* reason; // why no twin was found
  };
  bool FindTwin(uint32_t from_pkg,
                const GraphTile& tile,
                const DirectedEdge& edge,
                const GraphId& start_tile,
                uint32_t start_id,
                const midgard::PointLL& start,
                const midgard::PointLL& end,
                uint32_t to_pkg,
                Twin& twin) const;
  // Index of the node of package `to` that is the twin of node `node` (real id) of package
  // `from` at ll: mutually nearest within the join tolerance. -1 if none.
  int64_t TwinNode(uint32_t from, const GraphId& node, const midgard::PointLL& ll, uint32_t to) const;
  // Local index of the edge of package `to` leaving node `node` of the real tile `tile` along the
  // way, whose end node is the twin of `far` (real node id of package `from`, at far_ll) and whose
  // length is within the join tolerance of `length`. -1 if none, with the reason.
  int MatchEdge(uint32_t to,
                const GraphTile& tile,
                uint32_t node,
                uint64_t way_id,
                uint32_t length,
                uint32_t from,
                const GraphId& far,
                const midgard::PointLL& far_ll,
                const char*& reason) const;

  // Index of the package's node nearest to ll within the join tolerance in the given (real) tile,
  // or -1. Ties go to the lower index.
  int64_t NearestNode(uint32_t pkg, const GraphId& tile_base, const midgard::PointLL& ll) const;

  static constexpr int kUndecided = -2;
  int FindOwner(const midgard::PointLL& p) const;
  // Owner shared by every point of a box, or kUndecided.
  int BoxOwner(const midgard::AABB2<midgard::PointLL>& box) const;
  // Cached BoxOwner of a level-2 tile.
  int CachedTileOwner(int32_t tile_id) const;
  // Whether p lies in a level-2 tile that an ownership boundary may cross.
  bool NearOwnershipBoundary(const midgard::PointLL& p) const;
  // Package whose polygon is nearest to p (newest on ties); best/second receive the two smallest
  // distances in meters. -1 without packages.
  int NearestPackage(const midgard::PointLL& p, double& best, double& second) const;
  // Owner computed for one point.
  int PointOwner(const midgard::PointLL& p) const;

  // The package's own, unmodified tile with the given real id, or nullptr if it does not ship it.
  // Lives as long as the set.
  const GraphTile* RawTile(uint32_t pkg, const GraphId& real_base) const;
  // Whether the package owns the start node of one of its edges (real id).
  bool EdgeOwned(uint32_t pkg, const GraphId& edge) const;
  // Coordinates of a node in a package's own tiles.
  bool NodeLL(uint32_t pkg, const GraphId& node, midgard::PointLL& ll) const;
  // Logs, once per package and tile, that the package references a tile it does not ship.
  void WarnMissingTile(uint32_t pkg, const GraphId& tile) const;
  // Drops the pages of the package tars from this process's resident memory; they are read back
  // from the files when needed.
  void ReleasePages() const;

  std::vector<std::unique_ptr<Package>> packages_;
  // real tile base -> packages holding it, in slot (copy index) order
  std::unordered_map<uint64_t, std::vector<uint32_t>> copies_;
  // copy tile base -> (package, index of the tile in the package)
  std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>> by_copy_;
  // max distance in meters between a node and its twin in another package
  double join_tolerance_m_ = 0;
  uint64_t join_key_ = 0;
  bool overlays_ = false;
  // copy ids registered in the process, released with the set
  std::vector<CopyId> registered_;

  // level-2 tile id -> owner of the whole tile, or kUndecided; fine cell key -> owner of the whole
  // cell, or kUndecided. Both only ever grow, and hold values that depend only on the polygons.
  mutable std::shared_mutex owner_mutex_;
  mutable std::unordered_map<int32_t, int> tile_owner_;
  mutable std::unordered_map<uint64_t, int> cell_owner_;

  struct Counters {
    std::atomic<uint64_t> joined{0}, lost{0}, disabled{0}, clean_tiles{0}, id_rewrites{0},
        full_rewrites{0}, overlay_tiles{0};
  };
  mutable Counters counters_;
};

} // namespace baldr
} // namespace valhalla
