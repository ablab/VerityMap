//
// Created by Andrey Bzikadze on 03/31/22.
//

#include "kmer_index/indexed_contigs.hpp"

namespace veritymap::mapper {

class Mapper {
  const Config config;
  logging::Logger &logger;
  const size_t nthreads = 1;
  const matches::Matcher matcher;
  const scoring::DPScorer dp_scorer;
  const chaining::Chainer chainer;
  const RollingHash<Config::HashParams::htype> &hasher;

 private:
  [[nodiscard]] chaining::Chains MapSingleQueryStrand(const indexed_contigs::IndexedContigs &indexed_targets,
                                                      const Contig &query,
                                                      const dna_strand::Strand &query_strand) const {
    chaining::Chains chains;
    for (int i = 0; i < indexed_targets.Size(); ++i) {
      const matches::Matches matches = matcher.GetMatches(indexed_targets, i, query, query_strand);
      if (matches.size() < config.chaining_params.min_matches) {
        continue;
      }
      const auto [scores, backtracks] = dp_scorer.GetScores(matches);
      chaining::Chains new_chains =
          chainer.GetChains(indexed_targets.Contigs().at(i), query, query_strand, matches, scores, backtracks);
      for (chaining::Chain &chain : new_chains) { chains.emplace_back(std::move(chain)); }
    }
    return chains;
  }

  [[nodiscard]] std::optional<chaining::Chain> MapSingleQuery(
      const Contig &query, const indexed_contigs::IndexedContigs &indexed_targets) const {
    using score_type = typename Config::ChainingParams::score_type;

    chaining::Chains chains = MapSingleQueryStrand(indexed_targets, query, dna_strand::Strand::forward);
    chaining::Chains chains_r = MapSingleQueryStrand(indexed_targets, query, dna_strand::Strand::reverse);

    for (chaining::Chain &chain : chains_r) { chains.emplace_back(std::move(chain)); }

    if (chains.empty())
      return std::nullopt;

    int64_t top_range = chains.front().Range(config.common_params.k);
    if (top_range < config.chaining_params.min_chain_range)
      return std::nullopt;

    if (chains.size() == 1)
      return std::move(chains.front());

    auto pr_it = std::max_element(chains.begin(), chains.end(),
                                  [](const auto &lhs, const auto &rhs) { return lhs.score < rhs.score; });
    Config::ChainingParams::score_type sc_score = 0;
    for (auto it = chains.begin(); it != chains.end(); ++it) {
      if (it == pr_it)
        continue;
      if (it->score > sc_score)
        sc_score = it->score;
    }
    const double top_score_prop = static_cast<double>(sc_score) / pr_it->score;
    if (top_score_prop < config.chaining_params.max_top_score_prop)
      return std::move(*pr_it);
    return std::nullopt;
  }

 public:
  Mapper(const Mapper &mapper) = delete;
  Mapper(Mapper &&mapper) = delete;
  Mapper &operator=(const Mapper &mapper) = delete;
  Mapper &operator=(Mapper &&mapper) = delete;

  Mapper(const Config &config, logging::Logger &logger, const size_t nthreads,
         const RollingHash<Config::HashParams::htype> &hasher)
      : config{config},
        logger{logger},
        nthreads{nthreads},
        matcher{config.kmer_indexer_params, hasher},
        dp_scorer{config.common_params, config.chaining_params},
        chainer{config.common_params, config.chaining_params},
        hasher{hasher} {}

  void ParallelRun(const indexed_contigs::IndexedContigs &indexed_targets, const std::vector<Contig> &queries,
                   const std::filesystem::path &chains_fn, const std::filesystem::path &sam_fn,
                   const std::string &cmd) {
    std::mutex chainsMutex;

    std::ofstream chains_os(chains_fn);
    std::ofstream sam_os(sam_fn);
    for (const Contig &target : indexed_targets.Contigs()) {
      sam_os << "@SQ\tSN:" << target.id << "\tLN:" << target.seq.size() << "\n";
    }
    sam_os << "@PG\tID:VerityMap\tPN:VerityMap\tVN:2.0\tCL:" << cmd << "\n";

    std::function<void(const Contig &query)> align_read = [&indexed_targets, &chainsMutex, &chains_os, &sam_os,
                                                           this](const Contig &query) {
      std::optional<chaining::Chain> chain = MapSingleQuery(query, indexed_targets);

      if (chain.has_value()) {
        std::string sam_record =
            chaining::chain2samrecord(chain.value(), config.common_params, config.chain2sam_params);
        chainsMutex.lock();
        chains_os << chain.value();
        sam_os << sam_record << "\n";
        chainsMutex.unlock();
      }
    };

    process_in_parallel(queries, align_read, nthreads, true);
  }
};

}// namespace veritymap::mapper