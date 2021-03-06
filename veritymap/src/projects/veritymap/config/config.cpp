//
// Created by Andrey Bzikadze on 06/19/21.
//

#include "config.hpp"

#include <fstream>
#include <map>

using namespace veritymap;

Config Config::load_config_file(const std::filesystem::path& config_fn) {
  std::ifstream is(config_fn);
  std::string key, str_val;
  std::map<std::string, std::string> m;
  while (is >> key >> str_val) {
    VERIFY(not m.contains(key));
    m[key] = str_val;
  }
  using std::stoull, std::stoll, std::stod, std::stoi;
  Config::CommonParams common_params{stoull(m.at("k")), (bool) stoi(m.at("diploid"))};
  Config::HashParams hash_params{stoull(m.at("base"))};
  Config::KmerIndexerParams::ApproximateKmerIndexerParams aprx_kmer_indexer_params{
      stod(m.at("false_positive_probability")), stod(m.at("exp_base")), stoi(m.at("nhash")),
      stoull(m.at("chunk_size"))};
  Config::KmerIndexerParams::ApproximateCanonKmerIndexerParams aprx_canon_kmer_indexer_params{
      stod(m.at("false_positive_probability_canon")), stod(m.at("exp_base_canon")), stoi(m.at("nhash_canon")),
      stoull(m.at("chunk_size_canon"))};
  Config::KmerIndexerParams::ApproximateCanonSingleThreadKmerIndexerParams aprx_canon_single_thread_kmer_indexer_params{
      stod(m.at("false_positive_probability_canon_single_thread")), stod(m.at("exp_base_canon_single_thread")),
      stoi(m.at("nhash_canon_single_thread"))};
  Config::KmerIndexerParams kmer_indexer_params{
      stoull(m.at("min_uncovered_len")), stoull(m.at("max_rare_cnt_target")),
      // stoull(m.at("max_rare_cnt_query")),
      stoull(m.at("k_step_size")), stoull(m.at("k_window_size")), stod(m.at("window_regular_density")),
      Config::KmerIndexerParams::str2strategy(m.at("strategy")), aprx_kmer_indexer_params,
      aprx_canon_kmer_indexer_params, aprx_canon_single_thread_kmer_indexer_params,
      stod(m.at("careful_upper_bnd_cov_mult"))};
  Config::Chain2SAMParams::KSW2Params ksw_2_params{
      static_cast<int8_t>(stoi(m.at("match_score"))), static_cast<int8_t>(stoi(m.at("mis_score"))),
      static_cast<int8_t>(stoi(m.at("gapo"))), static_cast<int8_t>(stoi(m.at("gape")))};
  Config::ChainingParams chaining_params{stoull(m.at("min_matches")),
                                         stod(m.at("min_score")),
                                         stod(m.at("max_top_score_prop")),
                                         stoll(m.at("max_jump")),
                                         stod(m.at("misassembly_penalty_base")),
                                         stod(m.at("diff_penalty_mult")),
                                         stoull(m.at("min_chain_range")),
                                         stoull(m.at("max_supp_dist_diff")),
                                         stoi(m.at("min_uniq_kmers")),
                                         stod(m.at("match_score_unique")),
                                         stod(m.at("match_score_dup")),
                                         stod(m.at("match_score_rare"))};
  Config::Chain2SAMParams chain_2_sam_params{stod(m.at("min_end_ident")), ksw_2_params};

  return Config{common_params, hash_params, kmer_indexer_params, chaining_params, chain_2_sam_params};
}
