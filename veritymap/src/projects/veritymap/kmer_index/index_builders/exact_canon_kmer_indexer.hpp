//
// Created by Andrey Bzikadze on 06/02/22.
//

#pragma once

#include "kmer_index_builder.hpp"
#include "kmer_window.hpp"

namespace veritymap::kmer_index_builder::exact_canon {

class ExactCanonKmerIndexer : public AbstractKmerIndexBuilder {
 public:
  ExactCanonKmerIndexer(const RollingHash<Config::HashParams::htype> &hasher, const Config::CommonParams &common_params,
                        const Config::KmerIndexerParams &kmer_indexer_params, logging::Logger &logger)
      : AbstractKmerIndexBuilder{/*nthreads=*/1, hasher, common_params, kmer_indexer_params, logger} {}

  [[nodiscard]] kmer_index::KmerIndex Build(const std::vector<Contig> &contigs) const override {
    kmer_index::KmerIndex::KmerCounter full_counter, counter;

    for (const Contig &contig : contigs) {
      if (contig.size() < hasher.k) {
        continue;
      }
      KWH<Config::HashParams::htype> kwh(hasher, contig.seq, 0);
      while (true) {
        full_counter[kwh.hash()] += 1;
        if (!kwh.hasNext()) {
          break;
        }
        kwh = kwh.next();
      }
    }

    kmer_index::KmerIndex::Kmer2Pos kmer2pos;
    for (const Contig &contig : contigs) {
      auto &k2p{kmer2pos.emplace_back()};
      if (contig.size() < hasher.k + kmer_indexer_params.k_step_size) {
        continue;
      }
      KWH<Config::HashParams::htype> kwh(hasher, contig.seq, 0);
      int64_t latest_pos{0};
      for (kmer_index::kmer_window::KmerMinimizerWindow window(kwh, kmer_indexer_params.k_step_size, full_counter);
           kwh.hasNext();
           kwh = kwh.next(), window.Add(full_counter.at(kwh.hash()), kwh.hash(), kwh.get_fhash(), kwh.pos)) {
        const int64_t cur_pos = window.GetMinimizerPos();
        if (cur_pos != latest_pos) {
          auto [freq, hash, fhash] = window.GetMinimizer();
          if (freq <= kmer_indexer_params.max_rare_cnt_target) {
            counter[fhash] = full_counter.at(hash);
            k2p[fhash].emplace_back(kwh.pos);
            latest_pos = cur_pos;
          }
        }
      }
    }

    return {kmer2pos, counter, contigs};
  }
};

}// namespace veritymap::kmer_index_builder::exact_canon