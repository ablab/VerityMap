//
// Created by Andrey Bzikadze on 2/19/21.
//

#pragma once

#include <common/parallel.h>

#include <common/logging.hpp>
#include <mutex>
#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>

#include "chaining.hpp"
#include "config/config.hpp"
#include "dp_scoring.hpp"
#include "hash_utils.hpp"
#include "kmer_index/indexed_contigs.hpp"
// #include "kmer_index/kmer_index.hpp"
#include "kmer_index/target_indexer.hpp"
#include "mapper.hpp"
#include "matches.hpp"

namespace veritymap {

class VerityMap {
  const Config config;
  logging::Logger &logger;
  const bool only_index = false;
  const bool careful_mode = true;
  const size_t nthreads = 1;
  const RollingHash<Config::HashParams::htype> hasher;

 private:
 public:
  VerityMap(const Config &config, logging::Logger &logger, const bool only_index, const bool careful_mode,
            const size_t nthreads)
      : config{config},
        logger{logger},
        only_index{only_index},
        careful_mode{careful_mode},
        nthreads{nthreads},
        hasher{config.common_params.k, config.hash_params.base} {}

  VerityMap(const VerityMap &) = delete;
  VerityMap(VerityMap &&) = delete;
  VerityMap &operator=(const VerityMap &) = delete;
  VerityMap &operator=(VerityMap &&) = delete;

  void Map(const std::filesystem::path &target_path, const std::filesystem::path &queries_path,
           const std::filesystem::path &outdir, const std::string &cmd,
           const std::optional<std::filesystem::path> &index_path) {
    std::vector<Contig> targets{io::SeqReader(target_path).readAllContigs()};
    for (const Contig &target : targets) {
      logger.info() << "Target length " << target.seq.size() << ", name " << target.id << " from " << target_path
                    << std::endl;
    }

    std::vector<Contig> queries(io::SeqReader(queries_path).readAllContigs());
    logger.info() << "Queries from " << queries_path << ", total " << queries.size() << " sequences " << std::endl;

    kmer_index::TargetIndexer target_indexer(config.common_params, config.kmer_indexer_params, logger, hasher);
    const indexed_contigs::IndexedContigs indexed_targets =
        target_indexer.GetIndexedTargets(targets, queries, index_path, outdir, nthreads);

    {
      const auto no_solid_kmers_fn = outdir / "no_solid_kmers.bed";
      std::ofstream no_solid_kmers_os(no_solid_kmers_fn);
      indexed_targets.NoSolidRegions2Bed(config.kmer_indexer_params.min_uncovered_len, config.common_params.k,
                                         no_solid_kmers_os);
      logger.info() << "Finished exporting long (>= " << config.kmer_indexer_params.min_uncovered_len << " bp) "
                    << "regions without solid k-mers to " << no_solid_kmers_fn << std::endl;
    }

    if (only_index) {
      return;
    }

    const auto chains_fn = outdir / "chains.tsv";
    const auto sam_fn = outdir / "alignments.sam";

    logger.info() << "Computing chains and sam records..." << std::endl;
    mapper::Mapper mapper(config, logger, nthreads, hasher);
    mapper.ParallelRun(indexed_targets, queries, chains_fn, sam_fn, cmd);

    logger.info() << "Finished outputting chains to " << chains_fn << " and sam records to " << sam_fn << std::endl;
  }
};

}// End namespace veritymap