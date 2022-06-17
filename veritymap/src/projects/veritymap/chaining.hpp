//
// Created by Andrey Bzikadze on 2/22/21.
//

#pragma once

#include "cigar.hpp"
#include "config/config.hpp"
#include "ksw_align.hpp"
#include "matches.hpp"
#include "strand.hpp"

namespace veritymap::chaining {

struct Chain {
  const Contig &target;
  const Contig &query;
  const dna_strand::Strand query_strand;
  matches::Matches matches;
  Config::ChainingParams::score_type score;

  [[nodiscard]] size_t query_size() const { return query.size(); }

  [[nodiscard]] Config::ChainingParams::match_pos_type query_st() const { return matches.front().query_pos; }
  [[nodiscard]] Config::ChainingParams::match_pos_type query_en() const { return matches.back().query_pos; }
  [[nodiscard]] Config::ChainingParams::match_pos_type target_st() const { return matches.front().target_pos; }
  [[nodiscard]] Config::ChainingParams::match_pos_type target_en() const { return matches.back().target_pos; }

  Chain(const Contig &target, const Contig &query, const dna_strand::Strand query_strand, matches::Matches matches,
        const Config::ChainingParams::score_type score)
      : target{target},
        query{query},
        query_strand{query_strand},
        matches{std::move(matches)},
        score{score} {}

  Chain(const Chain &) = delete;
  Chain &operator=(Chain) = delete;
  Chain(Chain &&) = default;
  Chain &operator=(Chain &&) = delete;
};

inline bool operator<(const Chain &lhs, const Chain &rhs) { return lhs.score < rhs.score; }
inline bool operator>(const Chain &lhs, const Chain &rhs) { return operator<(rhs, lhs); }
inline bool operator<=(const Chain &lhs, const Chain &rhs) { return !operator>(lhs, rhs); }
inline bool operator>=(const Chain &lhs, const Chain &rhs) { return !operator<(lhs, rhs); }

std::ostream &operator<<(std::ostream &os, const Chain &chain) {
  const std::string strand = dna_strand::strand2str(chain.query_strand);
  os << strand << "Aln " << chain.query.id << " " << chain.target.id << " " << chain.query_st() << " "
     << chain.query_en()// TODO add k
     << " " << chain.query_size() << " " << chain.target_st() << " " << chain.target_en() << " score " << chain.score
     << " num kmers " << size(chain.matches) << "\n";
  os << "Chain\n" << chain.matches;
  return os;
}

using Chains = std::vector<Chain>;

std::ostream &operator<<(std::ostream &os, const Chains &chains) {
  for (const auto &chain : chains) { os << chain; }
  return os;
}

class Chainer {
  const Config::CommonParams common_params;
  const Config::ChainingParams chaining_params;

 public:
  Chainer(const Chainer &) = delete;
  Chainer(Chainer &&) = delete;
  Chainer &operator=(const Chainer &) = delete;
  Chainer &operator=(Chainer &&) = delete;

  Chainer(const Config::CommonParams &common_params, const Config::ChainingParams &chaining_params)
      : common_params{common_params},
        chaining_params{chaining_params} {}

  [[nodiscard]] Chains GetChains(const Contig &target, const Contig &query, const dna_strand::Strand &query_strand,
                                 const matches::Matches &matches,
                                 const std::vector<typename Config::ChainingParams::score_type> &scores,
                                 const std::vector<size_t> &backtracks) const {
    using score_type = typename Config::ChainingParams::score_type;

    std::vector<size_t> index;
    for (size_t i = 0; i < scores.size(); ++i) { index.emplace_back(i); }
    std::sort(index.begin(), index.end(), [&scores, &matches](const size_t &lhs, const size_t &rhs) {
      if (scores[lhs] != scores[rhs]) {
        return scores[lhs] > scores[rhs];
      }
      return matches[lhs].target_pos < matches[rhs].target_pos;
    });

    std::vector<int> end_chain(matches.size(), 0);// should not use std::vector<bool>
    size_t def_backtrack = std::numeric_limits<size_t>::max();
    for (unsigned long i : index) {
      if (backtracks[i] == def_backtrack) {
        end_chain[i] = 1;
      }
    }
    // for (size_t i = 0; i < scores.size(); ++i) {
    //   std::cout << i << " " << matches[i].target_pos << " " << matches[i].query_pos << " " << matches[i].target_freq << " "
    //             << scores[i] << " " << backtracks[i] << "\n";
    // }
    // std::cout << "\n";

    using ScoredMatches = std::pair<matches::Matches, score_type>;

    std::vector<ScoredMatches> scored_matches_vec;
    for (const size_t en : index) {
      if (end_chain[en]) {
        continue;
      }
      matches::Matches chain_matches;
      size_t st = en;
      while (not end_chain[st]) {
        chain_matches.push_back(matches[st]);
        end_chain[st] = 1;
        st = backtracks[st];
      }
      std::reverse(chain_matches.begin(), chain_matches.end());
      for (size_t i = 1; i < chain_matches.size(); ++i) {
        VERIFY(chain_matches[i].query_pos > chain_matches[i - 1].query_pos);
        VERIFY(chain_matches[i].target_pos > chain_matches[i - 1].target_pos);
      }
      VERIFY(not chain_matches.empty());
      const score_type score = scores[en] - scores[st];
      const size_t chain_range = [&chain_matches, this] {
        const size_t query_range = chain_matches.back().query_pos + common_params.k - chain_matches.front().query_pos;
        const size_t target_range =
            chain_matches.back().target_pos + common_params.k - chain_matches.front().target_pos;
        return std::max(query_range, target_range);
      }();
      const int uniq_kmers = [&chain_matches] {
        int uniq_kmers{0};
        for (const auto &match : chain_matches) {
          if (match.target_freq == 1) {
            ++uniq_kmers;
          }
        }
        return uniq_kmers;
      }();
      if ((score >= chaining_params.min_score) and (chain_range >= chaining_params.min_chain_range)
          and (uniq_kmers >= chaining_params.min_uniq_kmers)) {
        scored_matches_vec.emplace_back(std::move(chain_matches), score);
      }
    }

    std::sort(scored_matches_vec.begin(), scored_matches_vec.end(),
              [](const auto &lhs, const auto &rhs) { return lhs.second > rhs.second; });

    Chains chains_vec;
    for (auto &&scored_matches : scored_matches_vec) {
      chains_vec.emplace_back(target, query, query_strand, scored_matches.first, scored_matches.second);
    }

    return chains_vec;
  }
};

struct SemiInterval {
  // [ )
  long long t_st, t_en, q_st, q_en;
};
using SemiIntervals = std::vector<SemiInterval>;

SemiIntervals get_intervals(const Chain &chain, const size_t k) {
  const long long k_{static_cast<long long>(k)};
  SemiIntervals intervals;
  VERIFY(not chain.matches.empty());
  const auto &fst_match = chain.matches.front();

  const Sequence query_seq = chain.query_strand == dna_strand::Strand::forward ? chain.query.seq : chain.query.RC().seq;
  // if (fst_match.query_pos != 0) {
  intervals.push_back({std::max<long long>(0, fst_match.target_pos - fst_match.query_pos), fst_match.target_pos, 0,
                       fst_match.query_pos});
  // }
  for (auto it2 = chain.matches.cbegin(), it1 = it2++; it2 != chain.matches.cend(); ++it1, ++it2) {
    const auto &fst = *it1;
    const auto &snd = *it2;
    VERIFY((snd.query_pos >= fst.query_pos) and (snd.target_pos >= fst.target_pos));
    const size_t query_range = snd.query_pos - fst.query_pos;
    const size_t target_range = snd.target_pos - fst.target_pos;
    if ((query_range < k_) or (target_range < k_)) {
      VERIFY(query_range == target_range);
      continue;
    }
    if ((query_range == k_) and (target_range == k_)) {
      continue;
    }
    intervals.push_back({fst.target_pos + k_, snd.target_pos, fst.query_pos + k_, snd.query_pos});
  }
  const auto &lst_match = chain.matches.back();
  VERIFY(lst_match.query_pos + k <= chain.query_size());
  // if (lst_match.query_pos + k < chain.query_size()) {
  intervals.push_back(
      {lst_match.target_pos + k_,
       std::min<long long>(lst_match.target_pos + chain.query_size() - lst_match.query_pos, chain.target.size()),
       lst_match.query_pos + k_, static_cast<long long>(chain.query_size())});
  // }
  return intervals;
};

// TODO extract from here into cigar related source
std::vector<cigar_utils::Cigar> intervals2cigar(const SemiIntervals &intervals, const Chain &chain,
                                                const Sequence &query_seq, const size_t k, const int8_t match_score,
                                                const int8_t mis_score, const int8_t gapo, const int8_t gape,
                                                const double min_end_ident) {
  long long prev_pos = 0;
  std::vector<cigar_utils::Cigar> cigars;
  for (auto it = intervals.cbegin(); it != intervals.cend(); ++it) {
    const SemiInterval &interval = *it;
    const long long t_st = interval.t_st;
    const long long t_en = interval.t_en;
    const long long q_st = interval.q_st;
    const long long q_en = interval.q_en;

    if (q_st > prev_pos) {
      cigars.emplace_back(q_st - prev_pos, cigar_utils::CigarMode::M);
    }
    if ((t_st == t_en) and (q_st == q_en)) {
      // If the interval is empty it is either the very first one or the last one
      VERIFY((it == intervals.cbegin()) or (intervals.cend() - it == 1));
      cigars.emplace_back(0, cigar_utils::CigarMode::M);// CigarMode here does not matter
    } else if (t_st == t_en) {
      VERIFY(q_st != q_en);
      cigars.emplace_back(q_en - q_st, cigar_utils::CigarMode::I);
    } else if (q_st == q_en) {
      VERIFY(t_st != t_en);
      cigars.emplace_back(t_en - t_st, cigar_utils::CigarMode::D);
    } else {
      const Sequence target_substr = chain.target.seq.Subseq(t_st, t_en);
      const Sequence query_substr = query_seq.Subseq(q_st, q_en);
      cigar_utils::Cigar cigar_interval =
          ksw_align::align(target_substr, static_cast<int>(t_st), static_cast<int>(t_en), query_substr,
                           static_cast<int>(q_st), static_cast<int>(q_en), match_score, mis_score, gapo, gape);
      if ((it == intervals.cbegin()) or (intervals.cend() - it == 1)) {
        const double identity = cigar_interval.identity(target_substr, query_substr);
        // If low identity, prefer to soft clip
        if (identity <= min_end_ident) {
          cigars.emplace_back(q_en - q_st, cigar_utils::CigarMode::S);
        } else {
          cigars.emplace_back(std::move(cigar_interval));
        }
      } else {
        cigars.emplace_back(std::move(cigar_interval));
      }
    }
    prev_pos = q_en;
  }
  if (prev_pos + k <= chain.query_size()) {
    const Sequence target_end =
        chain.target.seq.Subseq(intervals.back().t_en, intervals.back().t_en + query_seq.size() - prev_pos);
    const Sequence query_suffix = query_seq.Subseq(prev_pos);
    VERIFY(target_end == query_suffix);
    cigars.emplace_back(query_suffix.size(), cigar_utils::CigarMode::M);
  }
  return cigars;
}

std::string chain2samrecord(const Chain &chain, const Config::CommonParams &common_params,
                            const Config::Chain2SAMParams &chain2sam_params) {
  const SemiIntervals intervals = get_intervals(chain, common_params.k);

  const Sequence query_seq = chain.query_strand == dna_strand::Strand::forward ? chain.query.seq : chain.query.RC().seq;

  const auto [cigar, left_trim] = [&intervals, &chain, &query_seq, &common_params, &chain2sam_params]() {
    cigar_utils::Cigar cigar;
    const Config::Chain2SAMParams::KSW2Params &ksw2p = chain2sam_params.ksw2_params;
    std::vector<cigar_utils::Cigar> cigars =
        intervals2cigar(intervals, chain, query_seq, common_params.k, ksw2p.match_score, ksw2p.mis_score, ksw2p.gapo,
                        ksw2p.gape, chain2sam_params.min_end_ident);
    for (cigar_utils::Cigar &cigar_ : cigars) { cigar.extend(std::move(cigar_)); }
    VERIFY(query_seq.size() == cigar.query_length());

    // If both front and end of cigars didn't consist of a single S then there will ==
    VERIFY(intervals.back().t_en - intervals.front().t_st >= cigar.target_length());

    const cigar_utils::CigarFragment &first_fragment = cigar.get_cigar_vec().front();
    size_t left_trim{0};
    if (first_fragment.mode != cigar_utils::CigarMode::S) {
      // If the first component is not soft clip, then need to remove deletions and soft clip the insertions
      const auto [left, right] = cigar.trim(cigar_utils::CigarMode::D);
      left_trim = left;
      cigar.soft_clip();
    } else {
      // If the first component is a soft clip, then need to adjust the starting target point
      left_trim = first_fragment.length;
    }
    return std::pair(cigar, left_trim);
  }();

  const Sequence target_seq =
      chain.target.seq.Subseq(intervals.front().t_st, intervals.front().t_st + cigar.target_length());
  const size_t start_pos = intervals.front().t_st + 1 + left_trim;
  const std::string read_flag = chain.query_strand == dna_strand::Strand::forward ? "0" : "16";
  std::stringstream s;
  s << chain.query.id << "\t" << read_flag << "\t" << chain.target.id << "\t" << start_pos << "\t60\t" << cigar
    << "\t*\t0\t0\t" << query_seq << "\t*\tNM:i:" << cigar.nmismatches(target_seq, query_seq);
  return s.str();
}

}// End namespace veritymap::chaining