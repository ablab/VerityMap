//
// Created by Andrey Bzikadze on 03/15/21.
//

#pragma once

#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>

#include "../config/config.hpp"
#include "../rolling_hash.hpp"
#include "common/logging.hpp"
#include "index_builders/exact_kmer_index_builder.hpp"
#include "indexed_contigs.hpp"
// #include "approx_canon_kmer_indexer.hpp"
// #include "approx_canon_kmer_indexer_single_thread.hpp"
// #include "approx_kmer_indexer.hpp"
// #include "bloom/bloom.hpp"
// #include "exact_canon_kmer_indexer.hpp"
// #include "include/sketch/ccm.h"

namespace veritymap::kmer_index {

using Counter = std::unordered_map<Config::HashParams::htype, size_t>;

class TargetIndexer {
  const Config::CommonParams common_params;
  const Config::KmerIndexerParams kmer_indexer_params;
  logging::Logger &logger;
  const RollingHash<Config::HashParams::htype> &hasher;

  KmerIndex ConstructIndex(const std::vector<Contig> &targets, const std::vector<Contig> &queries) {
    using namespace kmer_index_builder;
    if (kmer_indexer_params.strategy == Config::KmerIndexerParams::Strategy::exact) {
      logger.info() << "Getting exact kmer indexes..." << std::endl;
      ExactKmerIndexBuilder builder(hasher, common_params, kmer_indexer_params, logger);
      KmerIndex index = builder.Build(targets, queries);
      logger.info() << "Finished getting exact kmer indexes" << std::endl;
      return index;
    }
    /*
    if (kmer_indexer_params.strategy == Config::KmerIndexerParams::Strategy::approximate) {
      logger.info() << "Getting approximate kmer indexes..." << std::endl;

      const approx_kmer_indexer::ApproxKmerIndexer kmer_indexer(nthreads, hasher, common_params, kmer_indexer_params);
      KmerIndexes kmers_indexes = kmer_indexer.extract(targets, queries, logger);
      logger.info() << "Finished getting approximate kmer indexes" << std::endl;
      return kmers_indexes;
    }
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
                                                    const std::filesystem::path &outdir) {
    using htype = Config::HashParams::htype;

    const auto save_index_path = outdir / "kmer_indexes.tsv";
    auto index = [&index_path, &targets, &queries, this]() -> KmerIndex {
      if (index_path) {
        logger.info() << "Reading index from " << index_path.value() << "\n";
        return {targets, index_path.value()};
      }
      return ConstructIndex(targets, queries);
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