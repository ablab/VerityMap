//
// Created by Andrey Bzikadze on 06/14/22.
//

#pragma once

#include "../../config/config.hpp"
#include "../../rolling_hash.hpp"
#include "../kmer_index.hpp"

namespace veritymap::kmer_index_builder {

class AbstractKmerIndexBuilder {
 protected:
  size_t nthreads{1};
  const RollingHash<Config::HashParams::htype> &hasher;
  Config::CommonParams common_params;
  Config::KmerIndexerParams kmer_indexer_params;
  logging::Logger &logger;

 public:
  AbstractKmerIndexBuilder(const size_t nthreads, const RollingHash<Config::HashParams::htype> &hasher,
                           const Config::CommonParams &common_params,
                           const Config::KmerIndexerParams &kmer_indexer_params, logging::Logger &logger)
      : nthreads{nthreads},
        hasher{hasher},
        common_params{common_params},
        kmer_indexer_params{kmer_indexer_params},
        logger{logger} {}

  virtual kmer_index::KmerIndex Build(const std::vector<Contig> &queries, const std::vector<Contig> &target) const = 0;
};

};// End namespace veritymap::kmer_index_builder
