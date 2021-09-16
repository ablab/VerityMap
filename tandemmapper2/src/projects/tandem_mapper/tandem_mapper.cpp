//
// Created by Andrey Bzikadze on 2/19/21.
//

#include <iostream>
#include <iomanip>

#include <common/cl_parser.hpp>
#include <common/logging.hpp>

#include "version/version.hpp"
#include "tandem_mapper.hpp"
#include "config/config.hpp"
#include "config_dir_def/config_dir_def.hpp"


int main(int argc, char ** argv) {
    CLParser parser {{"output-dir=", "target=", "queries=", "threads=40",
                      "compress", "only-index", "index=none", "config=hifi"}, {},
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

    logging::LoggerStorage ls{output_dir, "tandem_mapper"};
    logging::Logger logger;
    const std::filesystem::path logfn = ls.newLoggerFile();
    logger.addLogFile(logfn);

    logger << "Log is written to " << logfn << std::endl;
    logger << "Git commit SHA1: " << tools::Version::GIT_SHA1 << std::endl;
    logger << "Git commit date: " << tools::Version::GIT_DATE << std::endl;

    auto time_point {std::chrono::system_clock::now()};
    std::time_t now = std::chrono::system_clock::to_time_t(time_point);
    logger << "Launch time: " << std::put_time(std::localtime(&now), "%c %Z") << std::endl;

    std::stringstream cmd_ss;
    for(size_t i = 0; i < argc; i++) {
        cmd_ss << argv[i] << " ";
    }
    const std::string cmd = cmd_ss.str();
    logger << "CMD: " << cmd << std::endl;

    const std::filesystem::path target_path =
            std::filesystem::canonical(parser.getValue("target"));
    const std::filesystem::path queries_path =
            std::filesystem::canonical(parser.getValue("queries"));

    bool to_compress = parser.getCheck("compress");
    bool only_index = parser.getCheck("only-index");


    const std::filesystem::path index_path = [&parser] {
        std::filesystem::path index_path = parser.getValue("index");
        if (index_path != "none") {
            index_path = std::filesystem::canonical(index_path);
        } else {
            index_path = "";
        }
        return index_path;
    }();


    const std::filesystem::path binary_path = argv[0];
    const std::filesystem::path config_fn = [&parser, &logger, &binary_path] {
        std::string config = parser.getValue("config");
        std::filesystem::path dirpath = binary_path.parent_path();
        if (config == "hifi") {
            return dirpath / "../../config/config_tm2_hifi.tsv";
        } else if (config == "ont") {
            return dirpath / "../../config/config_tm2_ont.tsv";
        }
        return static_cast<std::filesystem::path>(config);
    }();
    tandem_mapper::Config config = tandem_mapper::Config::load_config_file(config_fn);
    const auto config_out_fn = output_dir / "config.tsv";
    std::filesystem::copy_file(config_fn, config_out_fn,
                               std::filesystem::copy_options::overwrite_existing);
    logger.info() << "Config exported to " << config_out_fn << "\n";


    tandem_mapper::tandem_map(target_path, queries_path, output_dir,
                              to_compress, only_index, nthreads, logger, cmd, index_path, config);

    logger.info() << "Thank you for using TandemMapper2!" << std::endl;
}
