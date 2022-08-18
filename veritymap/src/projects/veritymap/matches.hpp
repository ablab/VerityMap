//
// Created by Andrey Bzikadze on 2/22/21.
//

#pragma once

#include "config/config.hpp"
#include "kmer_index/filter_rep_kmers.hpp"
#include "kmer_index/indexed_contigs.hpp"
#include "strand.hpp"

namespace veritymap::matches {

struct Match {
  Config::ChainingParams::match_pos_type target_pos{0};
  int32_t query_pos{0};
  uint8_t target_freq{0};// TODO Change to "is_unique"

  [[nodiscard]] bool is_unique() const { return target_freq == 1; }
  [[nodiscard]] bool is_dup() const { return target_freq == 2; }
};

inline bool operator<(const Match& lhs, const Match& rhs) { return lhs.target_pos < rhs.target_pos; }

inline std::ostream& operator<<(std::ostream& os, const Match& match) {
  os << match.query_pos << "\t" << match.target_pos << "\t" << static_cast<uint16_t>(match.target_freq) << "\n";
  return os;
}

using Matches = std::vector<Match>;

inline std::ostream& operator<<(std::ostream& os, const Matches& matches) {
  size_t prev_pos = 0;
  for (const auto& match : matches) {
    if (match.target_pos - prev_pos > 10) {// TODO: fix to k or k/2
      os << match;
      prev_pos = match.target_pos;
    }
  }
  return os;
}

class Matcher {
  const Config::KmerIndexerParams& config;
  const RollingHash<typename Config::HashParams::htype>& hasher;

 public:
  Matcher(const Matcher&) = delete;
  Matcher(Matcher&&) = delete;
  Matcher& operator=(const Matcher&) = delete;
  Matcher& operator=(Matcher&&) = delete;

  Matcher(const Config::KmerIndexerParams& config, const RollingHash<typename Config::HashParams::htype>& hasher)
      : config{config},
        hasher{hasher} {}

  [[nodiscard]] Matches GetMatches(const indexed_contigs::IndexedContigs& indexed_targets, const int64_t i,
                                   const Contig& query, const dna_strand::Strand& query_strand) const {
    Sequence seq = query_strand == dna_strand::Strand::forward ? query.seq : query.RC().seq;
    if (seq.size() < hasher.k) {
      return {};
    }

    // We are using approximate kmer detection
    // const double fpp{config.approximate_kmer_indexer_params.false_positive_probability};
    // BloomFilter rep_kmer_bf = veritymap::kmer_index::filter_rep_kmers::get_bloom_rep_kmers(seq, hasher, fpp);

    Matches matches;
    KWH<Config::HashParams::htype> kwh(hasher, seq, 0);
    const kmer_index::KmerIndex& index = indexed_targets.Index();
    while (true) {
      const Config::HashParams::htype hash = kwh.get_fhash();
      const std::vector<int64_t>* pos = index.GetPos(hash, i);
      const int64_t count_i = pos == nullptr ? 0 : pos->size();
      if (count_i) {
        const int64_t count64 = index.GetCount(hash);
        VERIFY(count64 >= count_i);
        VERIFY(count64 <= config.max_rare_cnt_target);
        VERIFY(count64 <= std::numeric_limits<uint8_t>::max());
        const auto count = static_cast<uint8_t>(count64);
        for (const int64_t tp : *pos) {
          VERIFY(kwh.pos < std::numeric_limits<int32_t>::max());
          matches.push_back(
              {static_cast<Config::ChainingParams::match_pos_type>(tp), static_cast<int32_t>(kwh.pos), count});
        }
      }
      if (not kwh.hasNext()) {
        break;
      }
      kwh = kwh.next();
    }
    std::sort(matches.begin(), matches.end());
    return matches;
  }
};

}// End namespace veritymap::matches