//
// Created by Andrey Bzikadze on 03/15/21.
//

#pragma once

#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>

#include "../config/config.hpp"
#include "../rolling_hash.hpp"
#include "common/logging.hpp"
#include "index_builders/approx_canon_kmer_indexer_single_thread.hpp"
#include "index_builders/approx_kmer_indexer_builder.hpp"
#include "index_builders/exact_canon_kmer_indexer.hpp"
#include "index_builders/exact_kmer_index_builder.hpp"
#include "indexed_contigs.hpp"
// #include "approx_canon_kmer_indexer_single_thread.hpp"
// #include "approx_kmer_indexer.hpp"
// #include "bloom/bloom.hpp"
// #include "include/sketch/ccm.h"

namespace veritymap::kmer_index {

using Counter = std::unordered_map<Config::HashParams::htype, size_t>;

class TargetIndexer {
  const Config::CommonParams common_params;
  const Config::KmerIndexerParams kmer_indexer_params;
  logging::Logger &logger;
  const RollingHash<Config::HashParams::htype> &hasher;

  KmerIndex ConstructIndex(const std::vector<Contig> &targets, const std::vector<Contig> &queries,
                           const int64_t nthreads) {
    using namespace kmer_index_builder;
    std::unique_ptr<AbstractKmerIndexBuilder> pbuilder;
    if (kmer_indexer_params.strategy == Config::KmerIndexerParams::Strategy::exact) {
      pbuilder = std::make_unique<exact::ExactKmerIndexBuilder>(hasher, common_params, kmer_indexer_params, logger);
    } else if (kmer_indexer_params.strategy == Config::KmerIndexerParams::Strategy::approximate) {
      pbuilder = std::make_unique<approx::ApproxKmerIndexBuilder>(nthreads, hasher, common_params, kmer_indexer_params,
                                                                  logger);
    } else if (kmer_indexer_params.strategy == Config::KmerIndexerParams::Strategy::exact_canon) {
      pbuilder =
          std::make_unique<exact_canon::ExactCanonKmerIndexer>(hasher, common_params, kmer_indexer_params, logger);
    } else if (kmer_indexer_params.strategy == Config::KmerIndexerParams::Strategy::approximate_canon) {
      pbuilder =
          std::make_unique<approx_canon::ApproxCanonKmerIndexer>(hasher, common_params, kmer_indexer_params, logger);
    }
    const std::string strategy_str = Config::KmerIndexerParams::strategy2str(kmer_indexer_params.strategy);
    logger.info() << "Getting kmer indexes. Strategy: " << strategy_str << std::endl;
    KmerIndex index = pbuilder->Build(targets);
    HighFreqUniqueKmersFilterer filterer(nthreads, hasher, common_params, kmer_indexer_params, logger);
    filterer.Filter(index, queries);
    logger.info() << "Finished getting kmer indexes" << std::endl;
    return index;
    /*
    if (kmer_indexer_params.strategy == Config::KmerIndexerParams::Strategy::approximate_canon) {
      logger.info() << "Getting approximate kmer indexes (canonical variant)..." << std::endl;
      const approx_canon_kmer_indexer::ApproxCanonKmerIndexer kmer_indexer(nthreads, hasher, common_params,
                                                                           kmer_indexer_params);
      KmerIndexes kmers_indexes = kmer_indexer.extract(targets, queries, logger);
      logger.info() << "Finished getting approximate (canonical variant) kmer indexes" << std::endl;
      return kmers_indexes;
    }
    if (kmer_indexer_params.strategy == Config::KmerIndexerParams::Strategy::exact_canon) {
      logger.info() << "Getting exact kmer indexes (canonical variant)..." << std::endl;
      const exact_canon_kmer_indexer::ExactCanonKmerIndexer kmer_indexer(hasher, common_params, kmer_indexer_params);
      KmerIndexes kmers_indexes = kmer_indexer.extract(targets, queries, logger);
      logger.info() << "Finished getting approximate (canonical variant) kmer indexes" << std::endl;
      return kmers_indexes;
    } else {
      VERIFY(kmer_indexer_params.strategy == Config::KmerIndexerParams::Strategy::approximate_canon_single_thread);
      const approx_canon_single_thread_kmer_indexer::ApproxCanonSingleThreadKmerIndexer kmer_indexer(
          hasher, common_params, kmer_indexer_params);
      KmerIndexes kmers_indexes = kmer_indexer.extract(targets, queries, logger);
      logger.info() << "Finished getting approximate (canonical variant, single thread) kmer indexes" << std::endl;
      return kmers_indexes;
    }
    */
  }

 public:
  TargetIndexer(const Config::CommonParams &common_params, const Config::KmerIndexerParams kmer_indexer_params,
                logging::Logger &logger, const RollingHash<Config::HashParams::htype> &hasher)
      : common_params{common_params},
        kmer_indexer_params{kmer_indexer_params},
        logger{logger},
        hasher{hasher} {}

  indexed_contigs::IndexedContigs GetIndexedTargets(const std::vector<Contig> &targets,
                                                    const std::vector<Contig> &queries,
                                                    const std::optional<std::filesystem::path> &index_path,
                                                    const std::filesystem::path &outdir, const int64_t nthreads) {
    using htype = Config::HashParams::htype;

    const auto save_index_path = outdir / "kmer_indexes.tsv";
    auto index = [&index_path, &targets, &queries, &nthreads, this]() -> KmerIndex {
      if (index_path) {
        logger.info() << "Reading index from " << index_path.value() << "\n";
        return {targets, index_path.value()};
      }
      return ConstructIndex(targets, queries, nthreads);
    }();

    std::ofstream index_os(save_index_path);
    index_os << index;
    index_os.close();
    logger.info() << "Kmer indexes are exported to " << save_index_path << std::endl;

    indexed_contigs::IndexedContigs indexed_targets(targets, hasher, index);
    indexed_targets.Summary(logger);

    return indexed_targets;
  }
};

}// End namespace veritymap::kmer_index