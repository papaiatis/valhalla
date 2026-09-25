// Writes the tile reference index ("tile_refs" in a mjolnir.packages entry) of a routing package:
// per tile, where its GraphIds lead. PackageSet uses it to serve tiles that need no rewrite
// straight from the memory-mapped package. The index depends only on the package, so it is built
// once, next to the package's tiles. See PackageSet::WriteTileRefs for its contents.
#include "argparse_utils.h"
#include "baldr/packageset.h"

#include <cxxopts.hpp>

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

using namespace valhalla;

int main(int argc, char** argv) {
  const auto program = std::filesystem::path(__FILE__).stem().string();
  std::string tile_extract, output;
  try {
    cxxopts::Options options(program, program + " " + VALHALLA_PRINT_VERSION +
                                          "\n\nBuilds the tile reference index of a package.\n");
    // clang-format off
    options.add_options()
      ("h,help", "Print this help message.")
      ("e,tile-extract", "Package tile tar.", cxxopts::value<std::string>(tile_extract))
      ("o,output", "Index file to write.", cxxopts::value<std::string>(output));
    // clang-format on
    auto result = options.parse(argc, argv);
    if (result.count("help") || tile_extract.empty() || output.empty()) {
      std::cout << options.help() << std::endl;
      return result.count("help") ? EXIT_SUCCESS : EXIT_FAILURE;
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  try {
    const auto references = baldr::PackageSet::WriteTileRefs(tile_extract, output);
    std::cout << "{\"references\": " << references << "}" << std::endl;
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
