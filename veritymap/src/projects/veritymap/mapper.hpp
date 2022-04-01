//
// Created by Andrey Bzikadze on 03/31/22.
//

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
  [[nodiscard]] chaining::Chains MapSingleQueryStrand(const kmer_index::IndexedContig &indexed_target,
                                                      const Contig &query,
                                                      const dna_strand::Strand &query_strand) const {
    const matches::Matches matches =
        matcher.GetMatches(indexed_target.get_contig(), indexed_target.get_kmer_index(), query, query_strand);
    if (matches.size() < config.chaining_params.min_matches) {
      return {};
    }

    const auto [scores, backtracks] = dp_scorer.GetScores(matches);

    chaining::Chains chains =
        chainer.GetChains(indexed_target.get_contig(), query, query_strand, matches, scores, backtracks);
    return chains;
  }

  [[nodiscard]] std::optional<chaining::Chain> MapSingleQuery(const Contig &query,
                                                              const kmer_index::IndexedContigs &indexed_targets) const {
    using score_type = typename Config::ChainingParams::score_type;

    chaining::Chains chains;
    for (const kmer_index::IndexedContig &indexed_target : indexed_targets) {
      chaining::Chains chains_f = MapSingleQueryStrand(indexed_target, query, dna_strand::Strand::forward);
      chaining::Chains chains_r = MapSingleQueryStrand(indexed_target, query, dna_strand::Strand::reverse);

      for (size_t i = 0; i < 2; ++i) {
        if (chains_f.size() > i) {
          chains.emplace_back(std::move(chains_f[i]));
        }
        if (chains_r.size() > i) {
          chains.emplace_back(std::move(chains_r[i]));
        }
      }
    }
    if (chains.empty())
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

  void ParallelRun(const kmer_index::IndexedContigs &indexed_targets, const std::vector<Contig> &queries,
                   const std::filesystem::path &chains_fn, const std::filesystem::path &sam_fn,
                   const std::string &cmd) {
    std::mutex chainsMutex;
    using TargetQuery = std::tuple<const kmer_index::IndexedContig *, const Contig *>;

    std::ofstream chains_os(chains_fn);
    std::ofstream sam_os(sam_fn);
    for (const kmer_index::IndexedContig &itarget : indexed_targets) {
      const Contig &target = itarget.get_contig();
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