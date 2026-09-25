// Diagnostics, not installed: loads every tile of a multi-package configuration (mjolnir.packages)
// and reports how the packages were joined: edges redirected across packages, joins lost, edges
// disabled, tiles served from overlays, and the time spent loading and rewriting tiles.
#include "argparse_utils.h"
#include "baldr/graphtile.h"
#include "baldr/packageset.h"

#include <boost/property_tree/ptree.hpp>
#include <cxxopts.hpp>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace valhalla;

int main(int argc, char** argv) {
  const auto program = std::filesystem::path(__FILE__).stem().string();
  boost::property_tree::ptree config;
  bool list_restrictions = false;
  try {
    cxxopts::Options options(program, program + " " + VALHALLA_PRINT_VERSION +
                                          "\n\nReports how independently built packages join.\n");
    // clang-format off
    options.add_options()
      ("h,help", "Print this help message.")
      ("v,version", "Print the version of this software.")
      ("c,config", "Path to the json configuration file.", cxxopts::value<std::string>())
      ("i,inline-config", "Inline JSON config", cxxopts::value<std::string>())
      ("r,restrictions", "List complex restrictions whose edges lie in more than one package.");
    // clang-format on
    auto result = options.parse(argc, argv);
    list_restrictions = result.count("restrictions") > 0;
    if (!parse_common_args(program, options, result, &config, "mjolnir.logging")) {
      return EXIT_SUCCESS;
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  const auto start = std::chrono::steady_clock::now();
  std::unique_ptr<baldr::PackageSet> packages;
  try {
    packages = baldr::PackageSet::FromConfig(config.get_child("mjolnir"));
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  if (!packages) {
    std::cerr << "No mjolnir.packages in config" << std::endl;
    return EXIT_FAILURE;
  }
  const auto opened = std::chrono::steady_clock::now();

  size_t tiles = 0, bytes = 0, mutable_bytes = 0;
  for (const auto& id : packages->AllTiles()) {
    if (auto tile = packages->LoadTile(id)) {
      ++tiles;
      const auto* header = tile->header();
      bytes += header->end_offset();
      // the parts the rewrite changes: header, directed edges, transitions, bins, restrictions
      mutable_bytes += sizeof(baldr::GraphTileHeader) +
                       header->directededgecount() * sizeof(baldr::DirectedEdge) +
                       header->transitioncount() * sizeof(baldr::NodeTransition) +
                       header->bin_offset(baldr::kBinCount - 1).second * sizeof(baldr::GraphId) +
                       (header->edgeinfo_offset() - header->complex_restriction_forward_offset());
    }
  }
  const auto loaded = std::chrono::steady_clock::now();

  if (list_restrictions) {
    // owner of an edge: the owner of its first shape point (shapes live in the edge's own tile)
    auto edge_owner = [&](const baldr::GraphId& edge, uint64_t& way) -> int {
      const auto tile = packages->LoadTile(edge.tile_base());
      if (!tile || edge.id() >= tile->header()->directededgecount()) {
        return -1;
      }
      const auto info = tile->edgeinfo(tile->directededge(edge));
      way = info.wayid();
      const auto shape = info.shape();
      return shape.empty() ? -1 : packages->Owner(shape.front());
    };
    size_t total = 0, crossing = 0;
    std::map<std::pair<int, size_t>, size_t> by_owner; // (first owner, number of owners) -> count
    for (const auto& id : packages->AllTiles()) {
      const auto tile = packages->LoadTile(id);
      if (!tile) {
        continue;
      }
      for (uint32_t i = 0; i < tile->header()->directededgecount(); ++i) {
        if (!tile->directededge(i)->end_restriction()) {
          continue;
        }
        const baldr::GraphId to(id.tileid(), id.level(), i);
        for (const auto& restriction : tile->GetComplexRestrictions(true, to, baldr::kAllAccess)) {
          ++total;
          std::vector<baldr::GraphId> edges{restriction.from_graphid()};
          restriction.WalkVias([&](const baldr::GraphId* via) {
            edges.push_back(*via);
            return baldr::WalkingVia::KeepWalking;
          });
          edges.push_back(restriction.to_graphid());
          std::set<int> owners;
          std::string ways;
          for (const auto& edge : edges) {
            uint64_t way = 0;
            owners.insert(edge_owner(edge, way));
            ways += " " + std::to_string(way);
          }
          ++by_owner[{*owners.begin(), owners.size()}];
          if (owners.size() > 1) {
            ++crossing;
            // one JSON line per crossing restriction, for routing tests: the shapes of its first
            // and last edge in travel direction
            auto shape_json = [&](const baldr::GraphId& edge) {
              const auto tile = packages->LoadTile(edge.tile_base());
              if (!tile || edge.id() >= tile->header()->directededgecount()) {
                return std::string("[]");
              }
              const auto* de = tile->directededge(edge);
              auto shape = tile->edgeinfo(de).shape();
              if (!de->forward()) {
                std::reverse(shape.begin(), shape.end());
              }
              std::string out = "[";
              for (size_t k = 0; k < shape.size(); ++k) {
                out += (k ? "," : "") + std::string("[") + std::to_string(shape[k].lng()) + "," +
                       std::to_string(shape[k].lat()) + "]";
              }
              return out + "]";
            };
            std::cout << "{\"restriction\": " << static_cast<int>(restriction.type())
                      << ", \"ways\": \"" << ways.substr(1)
                      << "\", \"from\": " << shape_json(restriction.from_graphid())
                      << ", \"to\": " << shape_json(restriction.to_graphid()) << "}" << std::endl;
            std::cerr << "Crossing restriction in tile "
                      << std::to_string(baldr::PackageSet::Real(id)) << " slot "
                      << baldr::PackageSet::CopyIndex(id) << " ways" << ways << " owners";
            for (auto o : owners) {
              std::cerr << " " << o;
            }
            std::cerr << std::endl;
          }
        }
      }
    }
    for (const auto& [key, count] : by_owner) {
      std::cerr << "owner " << key.first << " (owners " << key.second << "): " << count << std::endl;
    }
    std::cout << "{\"complex_restrictions\": " << total << ", \"crossing\": " << crossing << "}"
              << std::endl;
    return EXIT_SUCCESS;
  }

  const auto stats = packages->stats();
  const auto secs = [](auto a, auto b) { return std::chrono::duration<double>(b - a).count(); };
  std::cout << "{\"tiles\": " << tiles << ", \"tile_bytes\": " << bytes
            << ", \"mutable_bytes\": " << mutable_bytes << ", \"joined\": " << stats.joined
            << ", \"lost\": " << stats.lost << ", \"disabled\": " << stats.disabled
            << ", \"clean_tiles\": " << stats.clean_tiles
            << ", \"id_rewrites\": " << stats.id_rewrites
            << ", \"full_rewrites\": " << stats.full_rewrites
            << ", \"overlay_tiles\": " << stats.overlay_tiles
            << ", \"uses_overlays\": " << (packages->uses_overlays() ? "true" : "false")
            << ", \"open_seconds\": " << secs(start, opened)
            << ", \"load_rewrite_seconds\": " << secs(opened, loaded) << "}" << std::endl;
  return EXIT_SUCCESS;
}
