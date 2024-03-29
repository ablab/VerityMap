//
// Created by Andrey Bzikadze on 03/31/21.
//

#pragma once

#include <ctime>

#include "../../config/config.hpp"
#include "../../rolling_hash.hpp"
#include "common/coverage_utils.hpp"
#include "kmer_filter.hpp"
#include "kmer_index_builder.hpp"
#include "kmer_window.hpp"

namespace veritymap::kmer_index_builder::approx {

class ApproxKmerIndexBuilder : public AbstractKmerIndexBuilder {
  struct HashPosType {
    Config::HashParams::htype fhash{0};
    size_t pos{0};
    kmer_index::kmer_filter::KmerType kmer_type{kmer_index::kmer_filter::KmerType::banned};
  };

  std::vector<size_t> BinHashesInChunk(std::vector<std::vector<HashPosType>> &hashes_pos,
                                       const kmer_index::kmer_filter::KmerFilter &kmer_filter,
                                       KWH<Config::HashParams::htype> &kwh, const int64_t ctg_ind) const {
    std::vector<size_t> sizes(nthreads, 0);
    auto process_chunk = [this, &hashes_pos, &sizes, &kmer_filter, &ctg_ind](size_t i) {
      std::vector<HashPosType> &hashes_pos_th = hashes_pos[i];
      const size_t size = sizes[i];
      for (int j = 0; j < size; ++j) {
        const Config::HashParams::htype fhash = hashes_pos_th[j].fhash;
        hashes_pos_th[j].kmer_type =
            kmer_filter.GetKmerType(ctg_ind, fhash, i, kmer_indexer_params.max_rare_cnt_target);
      }
    };

    const size_t chunk_size = kmer_indexer_params.approximate_kmer_indexer_params.chunk_size;

    for (size_t cnt = 0; cnt < chunk_size; ++cnt) {
      const Config::HashParams::htype fhash = kwh.get_fhash();
      const Config::HashParams::htype rhash = kwh.get_rhash();
      const size_t ithread = ((fhash * rhash) % (2 * nthreads)) / 2;

      if (hashes_pos[ithread].size() == sizes[ithread]) {
        hashes_pos[ithread].push_back({fhash, kwh.pos});
      } else {
        hashes_pos[ithread][sizes[ithread]] = {fhash, kwh.pos};
      }
      ++sizes[ithread];
      if (not kwh.hasNext()) {
        break;
      }
      kwh = kwh.next();
    }
    std::vector<std::thread> threads(nthreads);
    for (size_t i = 0; i < threads.size(); ++i) { threads[i] = std::thread(process_chunk, i); }
    for (auto &thread : threads) { thread.join(); }

    return sizes;
  }

  void GetKmerIndex(const Contig &contig, const kmer_index::kmer_filter::KmerFilter &kmer_filter,
                    kmer_index::KmerIndex::KmerCounter &kmer_counter,
                    kmer_index::KmerIndex::Kmer2PosSingle &kmer2pos_single, const int64_t contig_ind) const {
    if (contig.size() < hasher.k) {
      return;
    }

    std::vector<std::vector<HashPosType>> hashes_pos(nthreads);

    KWH<Config::HashParams::htype> kwh({hasher, contig.seq, 0});
    const size_t window_size = kmer_indexer_params.k_window_size;
    const size_t step_size = kmer_indexer_params.k_step_size;
    std::vector<std::tuple<size_t, Config::HashParams::htype, bool>> pos_hash_uniq;
    kmer_index::kmer_window::KmerWindow kmer_window(window_size, pos_hash_uniq);
    while (true) {
      logger.info() << "Pos = " << kwh.pos << "\n";
      logger.info() << "Running jobs for chunk \n";
      std::vector<size_t> sizes = BinHashesInChunk(hashes_pos, kmer_filter, kwh, contig_ind);

      logger.info() << "Preparing kmer positions for sort \n";
      for (size_t ithread = 0; ithread < nthreads; ++ithread) {
        const auto &hashes_pos_th = hashes_pos[ithread];
        for (size_t j = 0; j < sizes[ithread]; ++j) {
          const auto &hash_pos = hashes_pos_th[j];
          if (hash_pos.kmer_type == kmer_index::kmer_filter::KmerType::unique
              or hash_pos.kmer_type == kmer_index::kmer_filter::KmerType::rare) {
            const bool is_unique = hash_pos.kmer_type == kmer_index::kmer_filter::KmerType::unique;
            pos_hash_uniq.emplace_back(hash_pos.pos, hash_pos.fhash, is_unique);
          }
        }
      }

      logger.info() << "Sorting kmer positions \n";
      std::sort(pos_hash_uniq.begin(), pos_hash_uniq.end());

      logger.info() << "Extending kmer index \n";
      kmer_window.Reset();
      auto it = pos_hash_uniq.begin();
      for (; it != pos_hash_uniq.end(); ++it) {
        const auto [pos, hash, is_unique] = *it;
        kmer_window.Inc();
        if (kwh.hasNext() and kwh.pos - pos < window_size / 2) {
          break;
        }
        if ((kmer_window.RegularFrac() < kmer_indexer_params.window_regular_density) or (pos % step_size == 0)) {
          ++kmer_counter[hash];
          kmer2pos_single[hash].push_back(pos);
        }
      }
      pos_hash_uniq = {it, pos_hash_uniq.end()};

      logger.info() << "Finished working with the chunk \n";
      if (not kwh.hasNext()) {
        break;
      }
    }
  }

  void FilterFalsePositives(kmer_index::KmerIndex::Kmer2Pos &kmer2pos,
                            kmer_index::KmerIndex::KmerCounter &kmer_counter) const {
    for (auto it = begin(kmer_counter); it != end(kmer_counter);) {
      if (it->second > kmer_indexer_params.max_rare_cnt_target) {
        for (kmer_index::KmerIndex::Kmer2PosSingle &kmer2pos_single : kmer2pos) { kmer2pos_single.erase(it->first); }
        it = kmer_counter.erase(it);// previously this was something like m_map.erase(it++);
      } else
        ++it;
    }
  }
  [[nodiscard]] kmer_index::KmerIndex GetKmerIndex(const std::vector<Contig> &contigs,
                                                   const kmer_index::kmer_filter::KmerFilter &kmer_filter) const {
    kmer_index::KmerIndex::KmerCounter kmer_counter;
    kmer_index::KmerIndex::Kmer2Pos kmer2pos;
    for (auto it = contigs.cbegin(); it != contigs.cend(); ++it) {
      const Contig &contig{*it};
      logger.info() << "Creating index for contig " << contig.id << "\n";
      kmer_index::KmerIndex::Kmer2PosSingle &kmer2pos_single = kmer2pos.emplace_back();
      GetKmerIndex(contig, kmer_filter, kmer_counter, kmer2pos_single, it - contigs.cbegin());
    }
    logger.info() << "Filtering potential false positive solid k-mers...\n";
    FilterFalsePositives(kmer2pos, kmer_counter);
    logger.info() << "Finished filtering\n";
    return {kmer2pos, kmer_counter, contigs};
  }

 public:
  ApproxKmerIndexBuilder(const int64_t nthreads, const RollingHash<Config::HashParams::htype> &hasher,
                         const Config::CommonParams &common_params,
                         const Config::KmerIndexerParams &kmer_indexer_params, logging::Logger &logger)
      : AbstractKmerIndexBuilder{nthreads, hasher, common_params, kmer_indexer_params, logger} {}

  [[nodiscard]] kmer_index::KmerIndex Build(const std::vector<Contig> &contigs) const override {
    const kmer_index::kmer_filter::KmerFilterBuilder kmer_filter_builder{nthreads, hasher, common_params,
                                                                         kmer_indexer_params};
    logger.info() << "Creating kmer filter\n";
    const kmer_index::kmer_filter::KmerFilter kmer_filter = kmer_filter_builder.GetKmerFilter(contigs, logger);
    logger.info() << "Finished creating kmer filter. Using it to build kmer indexes\n";
    kmer_index::KmerIndex kmer_index = GetKmerIndex(contigs, kmer_filter);

    return kmer_index;
  }
};

}// namespace veritymap::kmer_index_builder::approx