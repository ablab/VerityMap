//
// Created by Andrey Bzikadze on 06/14/22.
//

#pragma once

#include <unordered_set>

#include "../../rolling_hash.hpp"
#include "kmer_index_builder.hpp"

namespace veritymap::kmer_index_builder::exact {

class ExactKmerIndexBuilder : public AbstractKmerIndexBuilder {
  [[nodiscard]] std::vector<kmer_index::KmerIndex::KmerCounter> GetCounters(const std::vector<Contig> &ctgs) const {
    std::vector<kmer_index::KmerIndex::KmerCounter> counters;
    for (auto it = ctgs.cbegin(); it != ctgs.cend(); ++it) {
      const Contig &ctg = *it;
      std::unordered_map<Config::HashParams::htype, int64_t> &counter{counters.emplace_back()};
      if (ctg.size() < hasher.k) {
        continue;
      }
      KWH<Config::HashParams::htype> kwh(hasher, ctg.seq, 0);
      while (true) {
        counter[kwh.get_fhash()] += 1;
        if (!kwh.hasNext()) {
          break;
        }
        kwh = kwh.next();
      }
    }
    return counters;
  }

  [[nodiscard]] std::unordered_map<Config::HashParams::htype, std::unordered_set<size_t>> Hash2Seqs(
      const std::vector<kmer_index::KmerIndex::KmerCounter> &counters) const {
    std::unordered_map<Config::HashParams::htype, std::unordered_set<size_t>> hash2seqs;
    for (auto it = counters.cbegin(); it != counters.cend(); ++it) {
      for (const auto &[hash, cnt] : *it) { hash2seqs[hash].insert(it - counters.cbegin()); }
    }
    return hash2seqs;
  }

  [[nodiscard]] kmer_index::KmerIndex Counter2Index(
      const std::vector<Contig> &ctgs,
      const std::unordered_map<Config::HashParams::htype, std::unordered_set<size_t>> &hash2seqs,
      const std::vector<kmer_index::KmerIndex::KmerCounter> &counters) const {
    kmer_index::KmerIndex::KmerCounter kmer_counter;
    kmer_index::KmerIndex::Kmer2Pos kmer2pos;
    for (auto it = ctgs.cbegin(); it != ctgs.cend(); ++it) {
      auto &k2p{kmer2pos.emplace_back()};
      const Contig &ctg = *it;
      const auto &counter = counters.at(it - ctgs.cbegin());
      if (ctg.size() < hasher.k) {
        continue;
      }
      KWH<Config::HashParams::htype> kwh(hasher, ctg.seq, 0);
      while (true) {
        Config::HashParams::htype fhash{kwh.get_fhash()};
        Config::HashParams::htype rhash{kwh.get_rhash()};
        if ((hash2seqs.at(fhash).size() == 1) and (not hash2seqs.contains(rhash))
            and (counter.at(fhash) <= kmer_indexer_params.max_rare_cnt_target)) {
          k2p[fhash].emplace_back(kwh.pos);
          ++kmer_counter[fhash];
        }
        if (!kwh.hasNext()) {
          break;
        }
        kwh = kwh.next();
      }
    }
    return {kmer2pos, kmer_counter, ctgs};
  }

 public:
  ExactKmerIndexBuilder(const RollingHash<Config::HashParams::htype> &hasher, const Config::CommonParams &common_params,
                        const Config::KmerIndexerParams &kmer_indexer_params, logging::Logger &logger)
      : AbstractKmerIndexBuilder{/*nthreads=*/1, hasher, common_params, kmer_indexer_params, logger} {}

  [[nodiscard]] kmer_index::KmerIndex Build(const std::vector<Contig> &contigs) const override {
    std::vector<kmer_index::KmerIndex::KmerCounter> counters = GetCounters(contigs);
    std::unordered_map<Config::HashParams::htype, std::unordered_set<size_t>> hash2seqs = Hash2Seqs(counters);
    return Counter2Index(contigs, hash2seqs, counters);
  }
};

}// End namespace veritymap::kmer_index_builder::exact
