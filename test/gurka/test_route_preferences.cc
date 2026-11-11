#include "baldr/graphreader.h"
#include "gurka.h"
#include "sif/costfactory.h"
#include "test.h"

#include <gtest/gtest.h>

#include <functional>

#if !defined(VALHALLA_SOURCE_DIR)
#define VALHALLA_SOURCE_DIR
#endif

using namespace valhalla;
using namespace valhalla::baldr;

namespace {

const std::vector<std::string> kProfiles = {"auto", "motorcycle"};

using options_t = std::unordered_map<std::string, std::string>;

// Prefix every option with the costing options path of the profile
options_t opts(const std::string& costing, const options_t& options) {
  options_t result;
  for (const auto& [key, value] : options)
    result["/costing_options/" + costing + "/" + key] = value;
  return result;
}

valhalla::Api route(const gurka::map& map,
                    const std::vector<std::string>& waypoints,
                    const std::string& costing,
                    const options_t& options = {}) {
  return gurka::do_action(valhalla::Options::route, map, waypoints, costing, opts(costing, options));
}

float leg_time(const valhalla::Api& api) {
  return api.directions().routes(0).legs(0).summary().time();
}

// The cost of an edge of the map, after `change` was applied to a copy of the edge, with the given
// costing options.
float edge_cost(
    const gurka::map& map,
    const std::string& costing,
    const options_t& options,
    const std::string& from,
    const std::string& to,
    const std::function<void(DirectedEdge&)>& change = [](DirectedEdge&) {}) {
  const auto api = route(map, {from, to}, costing, options);
  const auto cost = sif::CostFactory{}.Create(api.options());
  GraphReader reader(map.config.get_child("mjolnir"));
  const auto [edge_id, edge] = gurka::findEdgeByNodes(reader, map.nodes, from, to);
  auto copy = *edge;
  change(copy);
  return cost->EdgeCost(&copy, edge_id, reader.GetGraphTile(edge_id)).cost;
}

} // namespace

/*********************************************************************************************/
// Two ways of the same length between A and D: the top one is paved (preferred by default), the
// bottom one is a compacted road.
class PavedFactor : public ::testing::Test {
protected:
  static gurka::map map;

  static void SetUpTestSuite() {
    const std::string ascii_map = R"(
      A----B
      |    |
      C----D
    )";
    const gurka::ways ways = {
        {"AB", {{"highway", "secondary"}}},
        {"BD", {{"highway", "secondary"}}},
        {"AC", {{"highway", "secondary"}, {"surface", "compacted"}}},
        {"CD", {{"highway", "secondary"}, {"surface", "compacted"}}},
    };
    const auto layout = gurka::detail::map_to_coordinates(ascii_map, 100);
    map = gurka::buildtiles(layout, ways, {}, {}, "test/data/gurka_pref_paved");
  }
};
gurka::map PavedFactor::map = {};

TEST_F(PavedFactor, MovesTheRouteToTheUnpavedWay) {
  for (const auto& costing : kProfiles) {
    gurka::assert::raw::expect_path(route(map, {"A", "D"}, costing), {"AB", "BD"});
    gurka::assert::raw::expect_path(route(map, {"A", "D"}, costing, {{"paved_factor", "10"}}),
                                    {"AC", "CD"});
  }
}

TEST_F(PavedFactor, OnlyPavedEdgesCostMore) {
  for (const auto& costing : kProfiles) {
    const float paved = edge_cost(map, costing, {}, "A", "B");
    EXPECT_NEAR(edge_cost(map, costing, {{"paved_factor", "3"}}, "A", "B"), paved * 3, paved * 0.01);
    // compacted is not paved: the cost does not change (auto: the fixed unpaved penalty goes away)
    if (costing == "motorcycle")
      EXPECT_FLOAT_EQ(edge_cost(map, costing, {{"paved_factor", "3"}}, "A", "C"),
                      edge_cost(map, costing, {}, "A", "C"));
    else
      EXPECT_LT(edge_cost(map, costing, {{"paved_factor", "3"}}, "A", "C"),
                edge_cost(map, costing, {}, "A", "C"));
  }
}

TEST_F(PavedFactor, OutOfRangeIsClampedWithAWarning) {
  const auto api = route(map, {"A", "D"}, "auto", {{"paved_factor", "50"}});
  ASSERT_EQ(api.info().warnings().size(), 1);
  EXPECT_EQ(api.info().warnings(0).code(), 500);
  gurka::assert::raw::expect_path(api, {"AB", "BD"});
}

/*********************************************************************************************/
// Top: secondary, shorter in time; bottom: residential.
class ClassWeights : public ::testing::Test {
protected:
  static gurka::map map;

  static void SetUpTestSuite() {
    const std::string ascii_map = R"(
      A----B
      |    |
      C----D
    )";
    const gurka::ways ways = {
        {"AB", {{"highway", "secondary"}}},
        {"BD", {{"highway", "secondary"}}},
        {"AC", {{"highway", "residential"}}},
        {"CD", {{"highway", "residential"}}},
    };
    const auto layout = gurka::detail::map_to_coordinates(ascii_map, 100);
    map = gurka::buildtiles(layout, ways, {}, {}, "test/data/gurka_pref_class");
  }
};
gurka::map ClassWeights::map = {};

TEST_F(ClassWeights, WeightMovesTheRoute) {
  for (const auto& costing : kProfiles) {
    gurka::assert::raw::expect_path(route(map, {"A", "D"}, costing), {"AB", "BD"});
    gurka::assert::raw::expect_path(route(map, {"A", "D"}, costing, {{"secondary_factor", "8"}}),
                                    {"AC", "CD"});
  }
}

TEST_F(ClassWeights, EveryFactorIsAtLeastOne) {
  for (const auto& costing : kProfiles) {
    const float stock_secondary = edge_cost(map, costing, {}, "A", "B");
    const float stock_residential = edge_cost(map, costing, {}, "A", "C");
    // A weight below 1 makes the other classes more expensive; the class itself never gets cheaper.
    const options_t low = {{"secondary_factor", "0.5"}};
    EXPECT_FLOAT_EQ(edge_cost(map, costing, low, "A", "B"), stock_secondary);
    EXPECT_NEAR(edge_cost(map, costing, low, "A", "C"), stock_residential * 2,
                stock_residential * 0.01);

    // Every one of the five classes: set on a copy of the edge, the cost is the stock cost times
    // weight / smallest weight, never below the stock cost.
    const std::vector<std::pair<std::string, baldr::RoadClass>> classes = {
        {"motorway_factor", baldr::RoadClass::kMotorway},
        {"trunk_factor", baldr::RoadClass::kTrunk},
        {"primary_factor", baldr::RoadClass::kPrimary},
        {"secondary_factor", baldr::RoadClass::kSecondary},
        {"tertiary_factor", baldr::RoadClass::kTertiary},
    };
    for (const auto& [option, road_class] : classes) {
      auto as_class = [&](DirectedEdge& e) { e.set_classification(road_class); };
      const float base = edge_cost(map, costing, {}, "A", "B", as_class);
      // weight 4 on the class: the class costs 4x; weight 0.5: the class costs the same as stock
      // (the other classes cost twice)
      EXPECT_NEAR(edge_cost(map, costing, {{option, "4"}}, "A", "B", as_class), base * 4, base * 0.01)
          << option;
      EXPECT_NEAR(edge_cost(map, costing, {{option, "0.5"}}, "A", "B", as_class), base, base * 0.01)
          << option;
      EXPECT_NEAR(edge_cost(map, costing, {{option, "0.5"}}, "A", "B",
                            [](DirectedEdge& e) {
                              e.set_classification(baldr::RoadClass::kResidential);
                            }),
                  edge_cost(map, costing, {}, "A", "B",
                            [](DirectedEdge& e) {
                              e.set_classification(baldr::RoadClass::kResidential);
                            }) *
                      2,
                  base * 0.02)
          << option;
    }
  }
}

TEST_F(ClassWeights, OutOfRangeIsClampedWithAWarning) {
  const auto api = route(map, {"A", "D"}, "motorcycle", {{"motorway_factor", "0.01"}});
  ASSERT_EQ(api.info().warnings().size(), 1);
  EXPECT_EQ(api.info().warnings(0).code(), 500);
}

/*********************************************************************************************/
// Top: a straight way (A, B, D), 2 km. Bottom: a tight zigzag way with high curvature.
class StraightFactor : public ::testing::Test {
protected:
  static gurka::map map;

  static void SetUpTestSuite() {
    const std::string ascii_map = R"(
      A---------B---------D

        1 3 5 7
         2 4 6 8
    )";
    const gurka::ways ways = {
        {"AB", {{"highway", "secondary"}}},
        {"BD", {{"highway", "secondary"}}},
        {"A12345678D", {{"highway", "secondary"}}},
    };
    const auto layout = gurka::detail::map_to_coordinates(ascii_map, 100);
    map = gurka::buildtiles(layout, ways, {}, {}, "test/data/gurka_pref_straight");
  }
};
gurka::map StraightFactor::map = {};

TEST_F(StraightFactor, StraightEdgeCostsMoreByTheFormula) {
  GraphReader reader(map.config.get_child("mjolnir"));
  ASSERT_EQ(std::get<1>(gurka::findEdgeByNodes(reader, map.nodes, "A", "B"))->curvature(), 0u);
  for (const auto& costing : kProfiles) {
    const float stock = edge_cost(map, costing, {}, "A", "B");
    const options_t straight = {{"straight_factor", "2"}};
    // curvature 0: factor 1 + straight_factor; curvature 6 or more: factor 1
    EXPECT_NEAR(edge_cost(map, costing, straight, "A", "B"), stock * 3, stock * 0.01);
    auto curvature = [](uint32_t c) { return [c](DirectedEdge& e) { e.set_curvature(c); }; };
    EXPECT_NEAR(edge_cost(map, costing, straight, "A", "B", curvature(6)), stock, stock * 0.01);
    EXPECT_NEAR(edge_cost(map, costing, straight, "A", "B", curvature(3)), stock * 2, stock * 0.01);
  }
}

TEST_F(StraightFactor, MovesTheRouteToTheCurvyWay) {
  GraphReader reader(map.config.get_child("mjolnir"));
  const auto zigzag = gurka::findEdge(reader, map.nodes, "A12345678D", "D");
  EXPECT_GE(std::get<1>(zigzag)->curvature(), 6u);
  for (const auto& costing : kProfiles) {
    gurka::assert::raw::expect_path(route(map, {"A", "D"}, costing), {"AB", "BD"});
    gurka::assert::raw::expect_path(route(map, {"A", "D"}, costing,
                                          {{"straight_factor", "12"}, {"skip_shortcuts", "1"}}),
                                    {"A12345678D"});
  }
}

// Density cannot be set from a map, so the formula is tested on copies of one edge.
TEST_F(StraightFactor, UrbanFactorFormula) {
  for (const auto& costing : kProfiles) {
    auto density = [](uint32_t d) { return [d](DirectedEdge& e) { e.set_density(d); }; };
    const options_t urban = {{"urban_factor", "9"}};
    // no extra cost up to density 6, then 1 + factor * (density - 6) / 9
    EXPECT_FLOAT_EQ(edge_cost(map, costing, urban, "A", "B", density(6)),
                    edge_cost(map, costing, {}, "A", "B", density(6)));
    EXPECT_FLOAT_EQ(edge_cost(map, costing, urban, "A", "B", density(0)),
                    edge_cost(map, costing, {}, "A", "B", density(0)));
    for (const auto& [d, factor] : std::vector<std::pair<uint32_t, float>>{{15, 10}, {12, 7}}) {
      const float stock = edge_cost(map, costing, {}, "A", "B", density(d));
      EXPECT_NEAR(edge_cost(map, costing, urban, "A", "B", density(d)), stock * factor, stock * 0.01);
    }
  }
}

/*********************************************************************************************/
// The top route has a signal; the bottom route is longer.
class SignalPenalty : public ::testing::Test {
protected:
  static gurka::map map;

  static void SetUpTestSuite() {
    const std::string ascii_map = R"(
      A--B--C--D
      |        |
      E--F--G--H
    )";
    const gurka::ways ways = {
        {"AB", {{"highway", "secondary"}}},   {"BC", {{"highway", "secondary"}}},
        {"CD", {{"highway", "secondary"}}},   {"AE", {{"highway", "secondary"}}},
        {"EFGH", {{"highway", "secondary"}}}, {"HD", {{"highway", "secondary"}}},
    };
    // a signal on a junction node (B) is on the top route
    const gurka::nodes nodes = {{"B", {{"highway", "traffic_signals"}}}};
    const auto layout = gurka::detail::map_to_coordinates(ascii_map, 100);
    map = gurka::buildtiles(layout, ways, nodes, {}, "test/data/gurka_pref_signal_node");
  }
};
gurka::map SignalPenalty::map = {};

TEST_F(SignalPenalty, NodeSignalMovesTheRoute) {
  for (const auto& costing : kProfiles) {
    gurka::assert::raw::expect_path(route(map, {"A", "D"}, costing), {"AB", "BC", "CD"});
    gurka::assert::raw::expect_path(route(map, {"A", "D"}, costing, {{"signal_penalty", "300"}}),
                                    {"AE", "EFGH", "HD"});
    // the same in the other direction (the reverse search sees the signal from the other side)
    gurka::assert::raw::expect_path(route(map, {"D", "A"}, costing, {{"signal_penalty", "300"}}),
                                    {"HD", "EFGH", "AE"});
  }
}

TEST_F(SignalPenalty, PenaltyIsCostNotTime) {
  for (const auto& costing : kProfiles) {
    // a small penalty does not change the route, and not the time shown
    const auto stock = route(map, {"A", "D"}, costing);
    const auto small = route(map, {"A", "D"}, costing, {{"signal_penalty", "1"}});
    gurka::assert::raw::expect_path(small, {"AB", "BC", "CD"});
    EXPECT_EQ(leg_time(stock), leg_time(small));
  }
}

TEST_F(SignalPenalty, SignalInTheMiddleOfAnEdge) {
  // the same map, with the signal in the middle of the top way
  const std::string ascii_map = R"(
      A--1--D
      |     |
      E-----F
  )";
  const gurka::ways ways = {
      {"A1D", {{"highway", "secondary"}}},
      {"AEFD", {{"highway", "secondary"}}},
  };
  // a signal in the middle of an edge needs a direction, here the way is travelled A to D
  const gurka::nodes nodes = {
      {"1", {{"highway", "traffic_signals"}, {"traffic_signals:direction", "forward"}}}};
  const auto layout = gurka::detail::map_to_coordinates(ascii_map, 100);
  const auto signal_map =
      gurka::buildtiles(layout, ways, nodes, {}, "test/data/gurka_pref_signal_edge");
  GraphReader reader(signal_map.config.get_child("mjolnir"));
  ASSERT_TRUE(std::get<1>(gurka::findEdgeByNodes(reader, layout, "A", "D"))->traffic_signal());

  for (const auto& costing : kProfiles) {
    gurka::assert::raw::expect_path(route(signal_map, {"A", "D"}, costing), {"A1D"});
    gurka::assert::raw::expect_path(route(signal_map, {"A", "D"}, costing,
                                          {{"signal_penalty", "300"}}),
                                    {"AEFD"});
    EXPECT_EQ(leg_time(route(signal_map, {"A", "D"}, costing)),
              leg_time(route(signal_map, {"A", "D"}, costing, {{"signal_penalty", "1"}})));
  }
}

/*********************************************************************************************/
// A roundabout between P and Q (clockwise A B C D E F), a longer way around it below. The way from
// P to the roundabout is long, so that in P to Q the bidirectional search meets before the roundabout
// and the reverse search charges the entry; in Q to P the forward search charges it.
class RoundaboutPenalty : public ::testing::Test {
protected:
  static gurka::map map;

  static void SetUpTestSuite() {
    const std::string ascii_map = R"(
                                    B-----C
                                    |     |
      P-----------------------------A     D-----Q
      |                             |     |     |
      |                             F-----E     |
      |                                         |
      |                                         |
      |                                         |
      |                                         |
      |                                         |
      |                                         |
      |                                         |
      |                                         |
      |                                         |
      G-----------------------------------------H
    )";
    const gurka::ways ways = {
        {"PA", {{"highway", "secondary"}}},
        {"DQ", {{"highway", "secondary"}}},
        {"AB", {{"highway", "secondary"}, {"junction", "roundabout"}}},
        {"BC", {{"highway", "secondary"}, {"junction", "roundabout"}}},
        {"CD", {{"highway", "secondary"}, {"junction", "roundabout"}}},
        {"DE", {{"highway", "secondary"}, {"junction", "roundabout"}}},
        {"EF", {{"highway", "secondary"}, {"junction", "roundabout"}}},
        {"FA", {{"highway", "secondary"}, {"junction", "roundabout"}}},
        {"PG", {{"highway", "secondary"}}},
        {"GH", {{"highway", "secondary"}}},
        {"HQ", {{"highway", "secondary"}}},
    };
    const auto layout = gurka::detail::map_to_coordinates(ascii_map, 100);
    map = gurka::buildtiles(layout, ways, {}, {}, "test/data/gurka_pref_roundabout");
  }
};
gurka::map RoundaboutPenalty::map = {};

TEST_F(RoundaboutPenalty, MovesTheRouteInBothSearchDirections) {
  for (const auto& costing : kProfiles) {
    // The roundabout is entered at the start of the route (P to Q) or at the end of the search
    // (Q to P): the forward and the reverse search must charge the same penalty.
    gurka::assert::raw::expect_path(route(map, {"P", "Q"}, costing), {"PA", "AB", "BC", "CD", "DQ"});
    gurka::assert::raw::expect_path(route(map, {"P", "Q"}, costing, {{"roundabout_penalty", "1000"}}),
                                    {"PG", "GH", "HQ"});
    gurka::assert::raw::expect_path(route(map, {"Q", "P"}, costing), {"DQ", "DE", "EF", "FA", "PA"});
    gurka::assert::raw::expect_path(route(map, {"Q", "P"}, costing, {{"roundabout_penalty", "1000"}}),
                                    {"HQ", "GH", "PG"});
  }
}

TEST_F(RoundaboutPenalty, PenaltyIsCostNotTime) {
  for (const auto& costing : kProfiles) {
    const auto stock = route(map, {"P", "Q"}, costing);
    const auto small = route(map, {"P", "Q"}, costing, {{"roundabout_penalty", "1"}});
    gurka::assert::raw::expect_path(small, {"PA", "AB", "BC", "CD", "DQ"});
    EXPECT_EQ(leg_time(stock), leg_time(small));
  }
}

/*********************************************************************************************/
// Two countries (Canada north of 49 N, United States south of it). X and Y are in Canada: the short
// way between them goes through the United States, the way around stays in Canada.
class BorderPenalty : public ::testing::Test {
protected:
  static gurka::map map;

  static void SetUpTestSuite() {
    const std::string ascii_map = R"(
      M-----------N
      |           |
      |           |
      |           |
      X           Y
      |           |
      |           |
      S-----------T
    )";
    const gurka::ways ways = {
        {"XM", {{"highway", "secondary"}}}, {"MN", {{"highway", "secondary"}}},
        {"NY", {{"highway", "secondary"}}}, {"XS", {{"highway", "secondary"}}},
        {"ST", {{"highway", "secondary"}}}, {"TY", {{"highway", "secondary"}}},
    };
    const auto layout =
        gurka::detail::map_to_coordinates(ascii_map, 100, {-111.962238354, 49.003392362});
    map = gurka::buildtiles(layout, ways, {}, {}, "test/data/gurka_pref_border",
                            {{"mjolnir.admin",
                              {VALHALLA_SOURCE_DIR "test/data/language_admin.sqlite"}}});
  }
};
gurka::map BorderPenalty::map = {};

TEST_F(BorderPenalty, MovesTheRouteToTheDomesticWay) {
  GraphReader reader(map.config.get_child("mjolnir"));
  ASSERT_TRUE(std::get<1>(gurka::findEdgeByNodes(reader, map.nodes, "X", "S"))->ctry_crossing());
  ASSERT_TRUE(std::get<1>(gurka::findEdgeByNodes(reader, map.nodes, "T", "Y"))->ctry_crossing());
  for (const auto& costing : kProfiles) {
    gurka::assert::raw::expect_path(route(map, {"X", "Y"}, costing), {"XS", "ST", "TY"});
    gurka::assert::raw::expect_path(route(map, {"X", "Y"}, costing, {{"border_penalty", "600"}}),
                                    {"XM", "MN", "NY"});
    gurka::assert::raw::expect_path(route(map, {"Y", "X"}, costing, {{"border_penalty", "600"}}),
                                    {"NY", "MN", "XM"});
  }
}

/*********************************************************************************************/
// skip_shortcuts: the primary road has a shortcut A to C.
TEST(SkipShortcuts, FlagIsSetAndTheRouteIsTheSame) {
  const std::string ascii_map = R"(
      A---B--C
  )";
  const gurka::ways ways = {{"AB", {{"highway", "primary"}, {"surface", "paved"}}},
                            {"BC", {{"highway", "primary"}, {"surface", "paved"}}}};
  const auto layout = gurka::detail::map_to_coordinates(ascii_map, 50);
  const auto map = gurka::buildtiles(layout, ways, {}, {}, "test/data/gurka_pref_shortcut");

  GraphReader reader(map.config.get_child("mjolnir"));
  ASSERT_TRUE(
      reader.GetShortcut(std::get<0>(gurka::findEdgeByNodes(reader, layout, "A", "B"))).is_valid());

  for (const auto& costing : kProfiles) {
    const auto stock = route(map, {"A", "C"}, costing);
    const auto skip = route(map, {"A", "C"}, costing, {{"skip_shortcuts", "1"}});
    EXPECT_FALSE(sif::CostFactory{}.Create(stock.options())->SkipShortcuts());
    EXPECT_TRUE(sif::CostFactory{}.Create(skip.options())->SkipShortcuts());
    gurka::assert::raw::expect_path(skip, {"AB", "BC"});
    EXPECT_EQ(leg_time(stock), leg_time(skip));
  }
}

/*********************************************************************************************/
// Motorcycle exclude_unpaved (same map as test_exclude_unpaved_smoothness.cc)
TEST(MotorcycleExcludeUnpaved, UnpavedRoadsAreAvoided) {
  const std::string ascii_map = R"(
   E----F----G----H----I----A----J----K----L
                                 |    |
                                 |    |
                                 |    |
                                 |    |
                                 M----N
  )";
  const gurka::ways ways = {
      {"HI", {{"highway", "residential"}, {"surface", "unpaved"}}},
      {"IA", {{"highway", "residential"}, {"surface", "unpaved"}}},
      {"IJ", {{"highway", "residential"}}},
      {"JK", {{"highway", "residential"}, {"surface", "gravel"}}},
      {"KL", {{"highway", "residential"}}},
      {"JM", {{"highway", "residential"}}},
      {"MN", {{"highway", "residential"}}},
      {"NK", {{"highway", "residential"}}},
  };
  const auto layout = gurka::detail::map_to_coordinates(ascii_map, 100);
  const auto map = gurka::buildtiles(layout, ways, {}, {}, "test/data/gurka_pref_motorcycle_unpaved");

  const auto allowed =
      route(map, {"I", "L"}, "motorcycle", {{"exclude_unpaved", "0"}, {"use_trails", "1"}});
  gurka::assert::raw::expect_path(allowed, {"IJ", "JK", "KL"});
  const auto excluded = route(map, {"I", "L"}, "motorcycle", {{"exclude_unpaved", "1"}});
  gurka::assert::raw::expect_path(excluded, {"IJ", "JM", "MN", "NK", "KL"});
  // and in the other direction (the reverse search)
  gurka::assert::raw::expect_path(route(map, {"L", "I"}, "motorcycle", {{"exclude_unpaved", "1"}}),
                                  {"KL", "NK", "MN", "JM", "IJ"});
}
