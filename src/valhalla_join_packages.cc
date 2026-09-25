// Joins the packages of a multi-package configuration (mjolnir.packages) ahead of time, as an app
// would at install time, and writes each package's overlay (<name>.joined) to the output
// directory. With mjolnir.package_joined pointing at that directory, the packages route without
// rewriting tiles. See PackageSet::WriteOverlays.
#include "argparse_utils.h"
#include "baldr/packageset.h"

#include <boost/property_tree/ptree.hpp>
#include <cxxopts.hpp>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

using namespace valhalla;

int main(int argc, char** argv) {
  const auto program = std::filesystem::path(__FILE__).stem().string();
  boost::property_tree::ptree config;
  std::string out_dir;
  try {
    cxxopts::Options options(
        program,
        program + " " + VALHALLA_PRINT_VERSION +
            "\n\nJoins independently built packages (mjolnir.packages) ahead of time: writes, per "
            "package,\nthe tiles the join changes to <output>/<name>.joined. Point "
            "mjolnir.package_joined at the\noutput directory to serve the packages from them; "
            "overlays that do not match the packages\nare ignored and the packages are joined at "
            "runtime.\n");
    // clang-format off
    options.add_options()
      ("h,help", "Print this help message.")
      ("v,version", "Print the version of this software.")
      ("c,config", "Path to the json configuration file.", cxxopts::value<std::string>())
      ("i,inline-config", "Inline JSON config", cxxopts::value<std::string>())
      ("o,output", "Directory for the overlays. Defaults to mjolnir.package_joined.",
       cxxopts::value<std::string>());
    // clang-format on
    auto result = options.parse(argc, argv);
    if (!parse_common_args(program, options, result, &config, "mjolnir.logging")) {
      return EXIT_SUCCESS;
    }
    out_dir = result.count("output") ? result["output"].as<std::string>()
                                     : config.get<std::string>("mjolnir.package_joined", "");
    if (out_dir.empty()) {
      throw std::runtime_error("--output or mjolnir.package_joined is required");
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  try {
    const auto start = std::chrono::steady_clock::now();
    const auto report = baldr::PackageSet::WriteOverlays(config.get_child("mjolnir"), out_dir);
    const std::chrono::duration<double> seconds = std::chrono::steady_clock::now() - start;
    std::cout << "{\"seconds\": " << seconds.count() << ", \"join_key\": \"" << std::hex << report.key
              << std::dec << "\", \"joined\": " << report.stats.joined
              << ", \"lost\": " << report.stats.lost
              << ", \"full_rewrites\": " << report.stats.full_rewrites
              << ", \"id_rewrites\": " << report.stats.id_rewrites
              << ", \"clean_tiles\": " << report.stats.clean_tiles << ", \"overlays\": {";
    for (size_t i = 0; i < report.overlays.size(); ++i) {
      const auto& overlay = report.overlays[i];
      std::cout << (i ? ", " : "") << "\"" << overlay.name << "\": {\"tiles\": " << overlay.tiles
                << ", \"bytes\": " << overlay.bytes << "}";
    }
    std::cout << "}}" << std::endl;
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
