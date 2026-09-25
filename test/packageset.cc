#include "baldr/graphid.h"
#include "baldr/graphtile.h"
#include "baldr/packageset.h"
#include "baldr/tilehierarchy.h"

#include <boost/property_tree/ptree.hpp>
#include <gtest/gtest.h>

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
