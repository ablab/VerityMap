//
// Created by Andrey Bzikadze on 2/22/21.
//

#pragma once

#include "config/config.hpp"

namespace veritymap::scoring {

using ScoresBacktracks = std::pair<std::vector<Config::ChainingParams::score_type>, std::vector<size_t>>;

class DPScorer {
  const Config::CommonParams &common_params;
  const Config::ChainingParams &chaining_params;
  std::vector<double> diff_pens;

 public:
  DPScorer(const DPScorer &) = delete;
  DPScorer(DPScorer &&) = delete;
  DPScorer &operator=(const DPScorer &) = delete;
  DPScorer &operator=(DPScorer &&) = delete;

  DPScorer(const Config::CommonParams &common_params, const Config::ChainingParams &chaining_params)
      : common_params{common_params},
        chaining_params{chaining_params} {
    for (int64_t i = 0; i < chaining_params.max_supp_dist_diff; ++i) { diff_pens.emplace_back(0); }
    for (uint64_t i = 1; i < chaining_params.max_jump; ++i) {
      diff_pens.emplace_back(std::sqrt(std::sqrt((double) i)) - 1);
    }
  }

  [[nodiscard]] ScoresBacktracks GetScores(const matches::Matches &matches) const {
    using score_type = Config::ChainingParams::score_type;
    using match_pos_type = Config::ChainingParams::match_pos_type;

    std::vector<score_type> scores;
    std::vector<size_t> backtracks;

    size_t def_backtrack = std::numeric_limits<size_t>::max();

    for (auto it = matches.cbegin(); it != matches.cend(); ++it) {
      const matches::Match &match{*it};

      const score_type freq_weight = match.is_unique() ? chaining_params.match_score_unique
          : match.is_dup()                             ? chaining_params.match_score_dup
                                                       : chaining_params.match_score_rare;

      score_type score{freq_weight};
      size_t backtrack{def_backtrack};
      VERIFY(scores.size() == it - matches.cbegin());

      for (auto [it2, sc_it] = std::pair{std::make_reverse_iterator(it), scores.crbegin()};
           (it2 != matches.crend()) and (sc_it != scores.crend()); ++it2, ++sc_it) {

        const matches::Match &prev_match{*it2};

        const match_pos_type target_jump = match.target_pos - prev_match.target_pos - common_params.k;
        if (target_jump >= chaining_params.max_jump) {
          break;
        }

        const match_pos_type query_jump = match.query_pos - prev_match.query_pos - common_params.k;

        if ((std::min(query_jump, target_jump) == -common_params.k) or// no jump on either query or target
            ((std::min(query_jump, target_jump) < 0) and (query_jump != target_jump)) or// non-equal overlap
            (query_jump >= chaining_params.max_jump)) {                                 // excessive jump on query
          continue;
        }

        const match_pos_type jump_penalty = std::min(std::abs(query_jump), std::abs(target_jump));
        const match_pos_type dist_diff = std::abs(std::abs(query_jump) - std::abs(target_jump));

        const score_type diff_penalty = diff_pens[dist_diff];

        const score_type overlap_penalty =
            std::min<score_type>(1, static_cast<score_type>(query_jump + common_params.k) / common_params.k);

        if (overlap_penalty < 1) {
          VERIFY(diff_penalty == 0);
        }
        const score_type cur_score = *sc_it + freq_weight * overlap_penalty
            - std::min(chaining_params.diff_penalty_mult * diff_penalty, chaining_params.misassembly_penalty);

        if (cur_score > score) {
          score = cur_score;
          VERIFY(scores.size() - 1 >= sc_it - scores.rbegin());
          backtrack = scores.size() - 1 - (sc_it - scores.rbegin());
          VERIFY(backtrack < scores.size());
        }
      }
      VERIFY(score >= freq_weight);
      scores.emplace_back(score);
      if (score == freq_weight) {
        VERIFY(backtrack == def_backtrack);
      }
      backtracks.emplace_back(backtrack);
    }
    VERIFY(scores.size() == matches.size());

    return {scores, backtracks};
  }
};

}// End namespace veritymap::scoring