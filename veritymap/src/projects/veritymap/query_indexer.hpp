//
// Created by Andrey Bzikadze on 03/15/21.
//

#pragma once

#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>

#include "../config/config.hpp"
#include "../rolling_hash.hpp"
#include "bloom/bloom.hpp"
#include "include/sketch/ccm.h"
#include "sketch_contigs.hpp"

namespace veritymap::kmer_index {

kmer_index::IndexedContigs
get_indexed_queries(const std::vector<Contig>& queries,
                    const std::filesystem::path& outdir,
                    const RollingHash<Config::HashParams::htype>& hasher,
                    const size_t nthreads,
                    logging::Logger& logger,
                    const std::filesystem::path& index_path,
                    const Config::CommonParams& common_params,
                    const Config::KmerIndexerParams& kmer_indexer_params) {
  using htype = Config::HashParams::htype;
  const std::vector<KmerIndex> kmers_indexes = [&queries, &hasher, &common_params, &kmer_indexer_params,
                                                &index_path, &logger] {
    if (kmer_indexer_params.strategy == Config::KmerIndexerParams::Strategy::exact) {
      logger.info() << "Getting exact kmer indexes..." << std::endl;
      std::vector<KmerIndex> kmers_indexes =
          get_rare_kmers(queries, hasher, kmer_indexer_params.max_rare_cnt_target);
      logger.info() << "Finished getting exact kmer indexes" << std::endl;
      return kmers_indexes;
    } else {
      VERIFY(kmer_indexer_params.strategy == Config::KmerIndexerParams::Strategy::approximate)
      logger.info() << "Getting approximate kmer indexes..." << std::endl;

      std::vector<KmerIndex> kmers_indexes =
          sketch_contigs::get_rare_kmers_approx<htype>(
              queries, hasher, common_params, kmer_indexer_params);
      logger.info() << "Finished getting approximate kmer indexes" << std::endl;
      return kmers_indexes;
    }
  }();

  IndexedContigs indexed_targets;
  for (auto it = kmers_indexes.begin(); it != kmers_indexes.end(); ++it) {
    const Contig& target = targets.at(it - kmers_indexes.begin());
    indexed_targets.emplace_back(target, hasher, kmer_indexer_params.max_rare_cnt_target, *it);
  }
  indexed_query.get_kmer_index().size() << std::endl;

  return indexed_query;
}
}// End namespace veritymap::kmer_index