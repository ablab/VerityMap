//
// Created by Andrey Bzikadze on 2/22/21.
//

#pragma once

#include <unordered_set>

#include "../config/config.hpp"
#include "../hash_utils.hpp"
#include "../rolling_hash.hpp"

namespace veritymap::kmer_index_ {

using KmerIndex = std::unordered_map<Config::HashParams::htype, std::vector<size_t>>;
using KmerIndexes = std::vector<KmerIndex>;

using Counter = std::unordered_map<Config::HashParams::htype, size_t>;

using Counters = std::vector<Counter>;

class IndexedContig {
  const Contig& contig;
  const RollingHash<Config::HashParams::htype>& hasher;
  size_t max_rare_cnt;
  KmerIndex kmer_index;

 public:
  IndexedContig(const Contig& contig_, const RollingHash<Config::HashParams::htype>& hasher_,
                const size_t max_rare_cnt_, KmerIndex kmer_index = {})
      : contig{contig_},
        hasher{hasher_},
        max_rare_cnt{max_rare_cnt_},
        kmer_index{std::move(kmer_index)} {}

  [[nodiscard]] const Contig& get_contig() const { return contig; }
  [[nodiscard]] size_t get_contig_size() const { return contig.size(); }
  const RollingHash<Config::HashParams::htype>& get_hasher() const { return hasher; }
  [[nodiscard]] size_t get_max_rare_cnt() const { return max_rare_cnt; }
  const KmerIndex& get_kmer_index() const { return kmer_index; }
};

using IndexedContigs = std::vector<IndexedContig>;

}// End namespace veritymap::kmer_index_