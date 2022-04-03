//
// Created by Andrey Bzikadze on 10/26/21.
//

#pragma once

#include <cmath>

#include "config/config.hpp"

namespace veritymap::cms_utils {
struct CMSParams {
  int nbits{1};
  int l2sz{1};
  int64_t nhash{1};

  CMSParams(const Config::CommonParams &common_params, const Config::KmerIndexerParams &kmer_indexer_params,
            const uint64_t tot_len, const uint nthreads)
      : nbits{(int) ceil(log2(kmer_indexer_params.max_rare_cnt_target))},
        l2sz{(int) ceil(
            log2(std::exp(kmer_indexer_params.approximate_kmer_indexer_params.exp_base) * tot_len / nthreads))},
        nhash{kmer_indexer_params.approximate_kmer_indexer_params.nhash} {}
};

}// namespace veritymap::cms_utils