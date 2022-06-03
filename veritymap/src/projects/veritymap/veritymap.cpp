//
// Created by Andrey Bzikadze on 2/19/21.
//

#include "veritymap.hpp"

#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <iomanip>
#include <iostream>

#include "config/config.hpp"
#include "version/version.hpp"

int main(int argc, char** argv) {
  CLParser parser{{"output-dir=", "target=", "queries=none", "threads=40", "only-index", "careful", "diploid",
                   "index=none", "config=hifi"},
                  {},
                  {"o=output-dir", "t=threads"}};
  parser.parseCL(argc, argv);
  if (!parser.check().empty()) {
    std::cerr << "Incorrect parameters" << std::endl;
    std::cerr << parser.check() << std::endl;
    return 1;
  }

  size_t nthreads = std::stoi(parser.getValue("threads"));
  if (nthreads == 0) {
    std::cerr << "# threads can not be set 0" << std::endl;
    return 1;
  }

  const std::filesystem::path output_dir{parser.getValue("output-dir")};
  ensure_dir_existance(output_dir);

  logging::LoggerStorage ls{output_dir, "veritymap"};
  logging::Logger logger;
  const std::filesystem::path logfn = ls.newLoggerFile();
  logger.addLogFile(logfn);

  logger << "Log is written to " << logfn << std::endl;
  logger << "Git commit SHA1: " << tools::Version::GIT_SHA1 << std::endl;
  logger << "Git commit date: " << tools::Version::GIT_DATE << std::endl;

  auto time_point{std::chrono::system_clock::now()};
  std::time_t now = std::chrono::system_clock::to_time_t(time_point);
  logger << "Launch time: " << std::put_time(std::localtime(&now), "%c %Z") << std::endl;

  std::stringstream cmd_ss;
  for (size_t i = 0; i < argc; i++) { cmd_ss << argv[i] << " "; }
  const std::string cmd = cmd_ss.str();
  logger << "CMD: " << cmd << std::endl;

  const std::filesystem::path target_path = std::filesystem::canonical(parser.getValue("target"));

  auto get_path_w_def = [&parser](const std::string& parameter) {
    std::filesystem::path path = parser.getValue(parameter);
    if (path != "none") {
      path = std::filesystem::canonical(path);
    } else {
      path = "";
    }
    return path;
  };
  const std::filesystem::path queries_path = get_path_w_def("queries");

  bool only_index = parser.getCheck("only-index");
  bool careful_mode = parser.getCheck("careful");
  if (careful_mode and queries_path == "") {
    std::cerr << "Cannot use careful mode if no queries are provided\n";
    return 1;
  }

  const std::filesystem::path index_path = get_path_w_def("index");

  const std::filesystem::path binary_path = argv[0];
  const std::filesystem::path config_fn = [&parser, &logger, &binary_path] {
    std::string config = parser.getValue("config");
    std::filesystem::path dirpath = binary_path.parent_path();
    if (config == "hifi") {
      return dirpath / "config/config_tm2_hifi.tsv";
    } else if (config == "ont") {
      return dirpath / "config/config_tm2_ont.tsv";
    }
    return static_cast<std::filesystem::path>(config);
  }();
  veritymap::Config config = veritymap::Config::load_config_file(config_fn);
  // bool diploid_mode = parser.getCheck("diploid");
  // if (diploid_mode) {
  //   // TODO refactor this out and modify config before copying into the output file
  //   config.common_params.diploid = true;
  //   config.kmer_indexer_params.strategy = veritymap::Config::KmerIndexerParams::Strategy::approximate_canon;
  // }

  const auto config_out_fn = output_dir / "config.tsv";
  std::filesystem::copy_file(config_fn, config_out_fn, std::filesystem::copy_options::overwrite_existing);
  logger.info() << "Config exported to " << config_out_fn << "\n";

  veritymap::VerityMap mapper(config, logger, only_index, careful_mode, nthreads);
  mapper.Map(target_path, queries_path, output_dir, cmd, index_path);

  logger.info() << "Thank you for using VerityMap!" << std::endl;
}
