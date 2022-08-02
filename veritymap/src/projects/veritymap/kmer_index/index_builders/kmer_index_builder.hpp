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
  int64_t nthreads{1};
  const RollingHash<Config::HashParams::htype> &hasher;
  Config::CommonParams common_params;
  Config::KmerIndexerParams kmer_indexer_params;
  logging::Logger &logger;

 public:
  AbstractKmerIndexBuilder(const int64_t nthreads, const RollingHash<Config::HashParams::htype> &hasher,
                           const Config::CommonParams &common_params,
                           const Config::KmerIndexerParams &kmer_indexer_params, logging::Logger &logger)
      : nthreads{nthreads},
        hasher{hasher},
        common_params{common_params},
        kmer_indexer_params{kmer_indexer_params},
        logger{logger} {}

  AbstractKmerIndexBuilder(const AbstractKmerIndexBuilder &) = delete;
  AbstractKmerIndexBuilder(AbstractKmerIndexBuilder &&) = delete;
  AbstractKmerIndexBuilder &operator=(const AbstractKmerIndexBuilder &) = delete;
  AbstractKmerIndexBuilder &operator=(AbstractKmerIndexBuilder &&) = delete;

  [[nodiscard]] virtual kmer_index::KmerIndex Build(const std::vector<Contig> &contigs) const = 0;
};

};// End namespace veritymap::kmer_index_builder
