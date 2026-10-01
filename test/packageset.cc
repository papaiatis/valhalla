#include "baldr/graphid.h"
#include "baldr/graphtile.h"
#include "baldr/packageset.h"
#include "baldr/tilehierarchy.h"

#include <boost/property_tree/ptree.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <limits>
#include <vector>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <type_traits>

using namespace valhalla::baldr;
using valhalla::midgard::PointLL;

namespace {

// Writes a .poly file into the test's temporary directory and returns its path.
std::string WritePoly(const std::string& name, const std::string& content) {
  const auto dir = std::filesystem::temp_directory_path() / "valhalla_packageset_test";
  std::filesystem::create_directories(dir);
  const auto path = (dir / (name + ".poly")).string();
  std::ofstream(path) << content;
  return path;
}

// A 10 x 10 degree square at the origin with a 2 x 2 degree hole in its middle.
const std::string kSquareWithHole = R"(square
outer
  0 0
  10 0
  10 10
  0 10
  0 0
END
!hole
  4 4
  6 4
  6 6
  4 6
END
END
)";

TEST(RegionPolygon, ContainsWithHole) {
  const RegionPolygon polygon(WritePoly("hole", kSquareWithHole));
  EXPECT_TRUE(polygon.Contains({1, 1}));
  EXPECT_TRUE(polygon.Contains({5, 1}));
  EXPECT_TRUE(polygon.Contains({5, 9}));
  EXPECT_FALSE(polygon.Contains({5, 5})) << "inside the hole";
  EXPECT_FALSE(polygon.Contains({11, 5}));
  EXPECT_FALSE(polygon.Contains({-1, 5}));
  EXPECT_FALSE(polygon.Contains({5, 10.5}));
}

TEST(RegionPolygon, MultipleOuterRings) {
  // two disjoint squares, the second one's ring left open (it is closed implicitly)
  const RegionPolygon polygon(WritePoly("two", R"(two
first
  0 0
  1 0
  1 1
  0 1
END
second
  5 0
  6 0
  6 1
  5 1
END
END
)"));
  EXPECT_TRUE(polygon.Contains({0.5, 0.5}));
  EXPECT_TRUE(polygon.Contains({5.5, 0.5}));
  EXPECT_FALSE(polygon.Contains({3, 0.5}));
  EXPECT_EQ(polygon.bbox().minx(), 0);
  EXPECT_EQ(polygon.bbox().maxx(), 6);
}

TEST(RegionPolygon, Distance) {
  const RegionPolygon polygon(WritePoly("hole_distance", kSquareWithHole));
  // 1 degree of latitude to the bottom edge
  EXPECT_NEAR(polygon.Distance({5, 1}), 110574.0, 1.0);
  // the hole's boundary counts: 1 degree of latitude below its lower edge
  EXPECT_NEAR(polygon.Distance({5, 3}), 110574.0, 1.0);
  // outside, 1 degree of longitude east of the right edge
  EXPECT_NEAR(polygon.Distance({11, 5}), 111320.0 * std::cos(5 * M_PI / 180), 1.0);
  // far away: beyond the cell search, every segment is checked
  EXPECT_NEAR(polygon.Distance({40, 5}), 30 * 111320.0 * std::cos(5 * M_PI / 180), 1.0);
}

// The distance to the nearest of the segments, in the local equirectangular projection around p
// that RegionPolygon::Distance uses, by looking at every segment.
double NearestSegmentMeters(const std::vector<std::array<double, 4>>& segments, double lng, double lat) {
  const double kx = 111320.0 * std::cos(lat * M_PI / 180);
  double best = std::numeric_limits<double>::max();
  for (const auto& [alng, alat, blng, blat] : segments) {
    const double ax = (alng - lng) * kx, ay = (alat - lat) * 110574.0;
    const double bx = (blng - lng) * kx, by = (blat - lat) * 110574.0;
    const double dx = bx - ax, dy = by - ay;
    const double len2 = dx * dx + dy * dy;
    const double t = len2 > 0 ? std::clamp(-(ax * dx + ay * dy) / len2, 0.0, 1.0) : 0.0;
    best = std::min(best, std::hypot(ax + t * dx, ay + t * dy));
  }
  return best;
}

// The search clips the rings of cells it visits to the cells that can hold a segment: it has to
// find the same nearest segment as looking at all of them, from inside, from the border, from just
// outside the polygon's cells on every side, and from far beyond the rings it searches.
TEST(RegionPolygon, DistanceEqualsTheNearestOfAllSegments) {
  const RegionPolygon polygon(WritePoly("hole_all_segments", kSquareWithHole));
  const std::vector<std::array<double, 4>> segments = {
      {0, 0, 10, 0}, {10, 0, 10, 10}, {10, 10, 0, 10}, {0, 10, 0, 0},
      {4, 4, 6, 4},  {6, 4, 6, 6},    {6, 6, 4, 6},    {4, 6, 4, 4},
  };
  size_t checked = 0;
  for (double lat : {-60.0, -20.0, -11.0, -0.3, 0.0, 2.2, 5.0, 7.9, 10.0, 10.4, 11.0, 20.0, 45.0, 80.0}) {
    for (double lng : {-120.0, -40.0, -11.0, -0.3, 0.0, 3.3, 5.0, 9.99, 10.0, 10.4, 12.0, 30.0, 100.0}) {
      const double expected = NearestSegmentMeters(segments, lng, lat);
      EXPECT_NEAR(polygon.Distance({lng, lat}), expected, 1e-6 * std::max(1.0, expected))
          << lng << "," << lat;
      ++checked;
    }
  }
  EXPECT_EQ(checked, 14u * 13u);
}

TEST(RegionPolygon, AntimeridianDoesNotThrow) {
  const RegionPolygon polygon(WritePoly("antimeridian", R"(antimeridian
ring
  179 -1
  -179 -1
  -179 1
  179 1
END
END
)"));
  // the ring is read as spanning the globe the long way round; queries must still work
  EXPECT_NO_THROW(polygon.Contains({180, 0}));
  EXPECT_NO_THROW(polygon.Distance({180, 0}));
}

TEST(RegionPolygon, MalformedFilesThrow) {
  EXPECT_THROW(RegionPolygon("/nonexistent/region.poly"), std::runtime_error);
  EXPECT_THROW(RegionPolygon(WritePoly("empty", "")), std::runtime_error);
  EXPECT_THROW(RegionPolygon(WritePoly("no_rings", "name\nEND\n")), std::runtime_error);
  EXPECT_THROW(RegionPolygon(WritePoly("unterminated", "name\nring\n 0 0\n 1 0\n 1 1\n")),
               std::runtime_error);
  EXPECT_THROW(RegionPolygon(WritePoly("bad_number", "name\nring\n 0 0\n 1 x\n 1 1\nEND\nEND\n")),
               std::runtime_error);
  EXPECT_THROW(RegionPolygon(WritePoly("out_of_range", "name\nring\n 0 0\n 1 0\n 1 91\nEND\nEND\n")),
               std::runtime_error);
  EXPECT_THROW(RegionPolygon(WritePoly("two_points", "name\nring\n 0 0\n 1 0\nEND\nEND\n")),
               std::runtime_error);
}

TEST(PackageSet, RealIdsMapToThemselves) {
  const GraphId real(1234, 2, 56);
  EXPECT_EQ(PackageSet::Real(real), real);
  EXPECT_EQ(PackageSet::CopyIndex(real), 0u);
  const GraphId invalid;
  EXPECT_FALSE(PackageSet::Real(invalid).is_valid());
  EXPECT_EQ(PackageSet::CopyIndex(invalid), 0u);
  // an id above the grid that no package set handed out is left alone
  const auto above = TileHierarchy::levels()[2].tiles.TileCount() + 1000;
  EXPECT_EQ(PackageSet::Real(GraphId(above, 2, 3)), GraphId(above, 2, 3));
  // transit tiles are not part of any package
  const GraphId transit(1234, 3, 0);
  EXPECT_EQ(PackageSet::Real(transit), transit);
}

TEST(PackageSet, SharedTilesAreThreadSafe) {
  boost::property_tree::ptree pt;
  boost::property_tree::ptree package;
  package.put("name", "region");
  package.put("tile_extract", "/nonexistent/region.tar");
  package.put("polygon", WritePoly("region", kSquareWithHole));
  boost::property_tree::ptree list;
  list.push_back({"", package});
  pt.put_child("packages", list);
#ifdef ENABLE_THREAD_SAFE_TILE_REF_COUNT
  // tiles of a package set are handed to every thread of the process
  EXPECT_TRUE((std::is_same_v<graph_tile_ptr, std::shared_ptr<const GraphTile>>));
  // the configuration is accepted and fails on the missing tar only
  try {
    PackageSet::FromConfig(pt);
    ADD_FAILURE() << "expected an error";
  } catch (const std::runtime_error& e) {
    EXPECT_NE(std::string(e.what()).find("/nonexistent/region.tar"), std::string::npos) << e.what();
  }
#else
  try {
    PackageSet::FromConfig(pt);
    ADD_FAILURE() << "expected an error";
  } catch (const std::runtime_error& e) {
    EXPECT_NE(std::string(e.what()).find("ENABLE_THREAD_SAFE_TILE_REF_COUNT"), std::string::npos)
        << e.what();
  }
#endif
}

TEST(PackageSet, NoPackagesConfigured) {
  boost::property_tree::ptree pt;
  EXPECT_EQ(PackageSet::FromConfig(pt), nullptr);
}

} // namespace
