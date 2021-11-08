//
// Created by Andrey Bzikadze on 06/19/21.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>

#include "sequences/verify.hpp"

namespace veritymap {
struct Config {
  struct CommonParams {
    size_t k;
  };
  CommonParams common_params;

  struct HashParams {
    using htype = uint64_t;
    size_t base;
  };
  HashParams hash_params;

  struct KmerIndexerParams {
    size_t min_uncovered_len;

    size_t max_rare_cnt_target;
    // size_t max_rare_cnt_query;
    // static_assert(max_rare_cnt_target <= std::numeric_limits<uint8_t>::max(),
    //               "Match stores frequency as uint8_t to economize memory usage");
    size_t k_step_size;
    size_t k_window_size;
    double window_unique_density;

    enum class Strategy { exact,
                          approximate };
    Strategy strategy;
    static std::string strategy2str(const Strategy& strategy) {
      return strategy == Strategy::exact ? "exact" : "approximate";
    }
    static Strategy str2strategy(const std::string& str) {
      VERIFY(str == "exact" or str == "approximate");
      if (str == "exact")
        return Strategy::exact;
      return Strategy::approximate;
    }

    struct ApproximateKmerIndexerParams {
      double false_positive_probability;
      double exp_base;
      int nhash;
      size_t chunk_size;
    };
    ApproximateKmerIndexerParams approximate_kmer_indexer_params;

    double careful_upper_bnd_cov_mult;
  };
  KmerIndexerParams kmer_indexer_params;

  struct ChainingParams {
    using match_pos_type = int64_t;
    using score_type = double;
    size_t min_matches;
    score_type min_score;
    double max_top_score_prop;
    match_pos_type max_jump;
    score_type misassembly_penalty;
    score_type diff_penalty_mult;
    size_t min_chain_range;
    size_t max_supp_dist_diff;
    // static_assert(min_chain_range / KmerIndexerParams::k_step_size >= min_matches);
    int min_uniq_kmers;

    score_type match_score_unique;
    score_type match_score_rare;
  };
  ChainingParams chaining_params;

  struct Chain2SAMParams {
    double min_end_ident;
    struct KSW2Params {
      int8_t match_score;
      int8_t mis_score;
      int8_t gapo;
      int8_t gape;
    };
    KSW2Params ksw2_params;
  };
  Chain2SAMParams chain2sam_params;

  static Config load_config_file(const std::filesystem::path& config_fn);
};

}// namespace veritymap