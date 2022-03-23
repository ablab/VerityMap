//
// Created by Andrey Bzikadze on 03/31/21.
//

#pragma once
#include <omp.h>

#include <ctime>

#include "../config/config.hpp"
#include "../rolling_hash.hpp"
#include "common/coverage_utils.hpp"
#include "kmer_filter.hpp"
#include "kmer_index.hpp"
#include "kmer_window.hpp"

namespace veritymap::kmer_index::approx_kmer_indexer {

template<typename htype>
class ApproxKmerIndexer {
  size_t nthreads{0};
  const RollingHash<htype> &hasher;
  Config::CommonParams common_params;
  Config::KmerIndexerParams kmer_indexer_params;

  struct HashPosType {
    htype fhash{0};
    size_t pos{0};
    kmer_filter::KmerType kmer_type{kmer_filter::KmerType::banned};
  };

  std::vector<size_t> BinHashesInChunk(std::vector<std::vector<HashPosType>> &hashes_pos,
                                       const kmer_filter::KmerFilter &kmer_filter,
                                       KWH<htype> &kwh,
                                       const size_t ctg_ind) const {
    std::vector<size_t> sizes(nthreads, 0);
    auto process_chunk = [this, &hashes_pos, &sizes, &kmer_filter, &ctg_ind](size_t i) {
      std::vector<HashPosType> &hashes_pos_th = hashes_pos[i];
      const size_t size = sizes[i];
      for (int j = 0; j < size; ++j) {
        const htype fhash = hashes_pos_th[j].fhash;
        hashes_pos_th[j].kmer_type = kmer_filter.GetKmerType(ctg_ind, fhash, i, kmer_indexer_params.max_rare_cnt_target);
      }
    };

    const size_t chunk_size = kmer_indexer_params.approximate_kmer_indexer_params.chunk_size;

    for (size_t cnt = 0; cnt < chunk_size; ++cnt) {
      const htype fhash = kwh.get_fhash();
      const htype rhash = kwh.get_rhash();
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
    for (size_t i = 0; i < threads.size(); ++i) {
      threads[i] = std::thread(process_chunk, i);
    }
    for (auto &thread : threads) {
      thread.join();
    }

    return sizes;
  }

  [[nodiscard]] KmerIndex GetKmerIndex(const Contig &contig,
                                       const kmer_filter::KmerFilter &kmer_filter,
                                       const size_t ctg_ind,
                                       logging::Logger &logger) const {
    if (contig.size() < hasher.k) {
      return {};
    }

    std::vector<std::vector<HashPosType>> hashes_pos(nthreads);

    KmerIndex kmer_index;
    KWH<htype> kwh({hasher, contig.seq, 0});
    const size_t window_size = kmer_indexer_params.k_window_size;
    kmer_window::KmerWindow kmer_window(window_size);
    const size_t step_size = kmer_indexer_params.k_step_size;
    while (true) {
      logger.info() << "Pos = " << kwh.pos << "\n";
      logger.info() << "Running jobs for chunk \n";
      std::vector<size_t> sizes = BinHashesInChunk(hashes_pos, kmer_filter, kwh, ctg_ind);

      logger.info() << "Preparing kmer positions for sort \n";
      std::vector<std::tuple<size_t, htype, bool>> pos_hash_uniq;
      for (size_t ithread = 0; ithread < nthreads; ++ithread) {
        const auto &hashes_pos_th = hashes_pos[ithread];
        for (size_t j = 0; j < sizes[ithread]; ++j) {
          const auto &hash_pos = hashes_pos_th[j];
          if (hash_pos.kmer_type == kmer_filter::KmerType::unique or hash_pos.kmer_type == kmer_filter::KmerType::rare) {
            const bool is_unique = hash_pos.kmer_type == kmer_filter::KmerType::unique;
            pos_hash_uniq.emplace_back(hash_pos.pos, hash_pos.fhash, is_unique);
          }
        }
      }

      logger.info() << "Sorting kmer positions \n";
      std::sort(pos_hash_uniq.begin(), pos_hash_uniq.end());

      logger.info() << "Extending kmer index \n";
      for (const auto &[pos, hash, is_unique] : pos_hash_uniq) {
        kmer_window.add(pos, is_unique);
        if ((kmer_window.unique_frac() < kmer_indexer_params.window_unique_density) or (pos % step_size == 0)) {
          kmer_index[hash].emplace_back(pos);
        }
      }

      logger.info() << "Finished working with the chunk \n";
      if (not kwh.hasNext()) {
        break;
      }
    }
    return kmer_index;
  }

  [[nodiscard]] KmerIndexes GetKmerIndexes(const std::vector<Contig> &contigs,
                                           const kmer_filter::KmerFilter &kmer_filter,
                                           logging::Logger &logger) const {
      KmerIndexes kmer_indexes;
      for (auto it = contigs.cbegin(); it!=contigs.cend(); ++it) {
          const Contig &contig{*it};
          logger.info() << "Creating index for contig " << contig.id << "\n";
          kmer_indexes.emplace_back(GetKmerIndex(contig,
                                                 kmer_filter,
                                                 it - contigs.cbegin(),
                                                 logger));
      }
      return kmer_indexes;
  }

    void BanHighFreqUniqueKmers(const std::vector<Contig> &contigs,
                                const std::vector<Contig> &readset,
                                KmerIndexes &kmer_indexes,
                                logging::Logger &logger) const {

        // ban unique k-mers in assembly that have unusually high coverage
        const double coverage
            {tools::common::coverage_utils::get_coverage(contigs, readset)};
        const uint max_read_freq = std::max(1.,
                                            ceil(kmer_indexer_params
                                                     .careful_upper_bnd_cov_mult
                                                     *coverage));

        Counter kmer_cnt;
        for (auto it = readset.begin(); it!=readset.end(); ++it) {
            logger.trace() << it - readset.begin() << " " << readset.size()
                           << "\n";
            const Contig &contig = *it;
            if (contig.size() < hasher.k) {
                continue;
            }
            KWH<htype> kwh(hasher, contig.seq, 0);
            while (true) {
                if (!kwh.hasNext()) {
                    break;
                }
                kwh = kwh.next();
                const htype fhash = kwh.get_fhash();
                const htype rhash = kwh.get_rhash();
                for (const htype hash : std::vector<htype>{fhash, rhash}) {
                    bool is_unique = false;
                    for (const KmerIndex &index : kmer_indexes) {
                        auto it = index.find(hash);
                        if (it!=index.end() and it->second.size()==1) {
                            is_unique = true;
                            break;
                        }
                    }
                    if (is_unique) {
                        kmer_cnt[hash] += 1;
                    }
                }
            }
        }

        uint64_t n{0};
        for (auto &[hash, cnt] : kmer_cnt) {
            if (cnt > max_read_freq) {
                for (KmerIndex &index : kmer_indexes) {
                    auto it = index.find(hash);
                    if (it!=index.end()) {
                        index.erase(it);
                        break;
                    }
                }
                ++n;
            }
        }
        logger.info() << "Filtered " << n << " high multiplicity k-mers\n";
    }

 public:
    ApproxKmerIndexer(const size_t nthreads,
                      const RollingHash<htype> &hasher,
                      const Config::CommonParams &common_params,
                      const Config::KmerIndexerParams &kmer_indexer_params)
        : nthreads{nthreads},
          hasher{hasher},
          common_params{common_params},
          kmer_indexer_params{
              kmer_indexer_params} {}

  ApproxKmerIndexer(const ApproxKmerIndexer &) = delete;
  ApproxKmerIndexer(ApproxKmerIndexer &&) = delete;
  ApproxKmerIndexer &operator=(const ApproxKmerIndexer &) = delete;
  ApproxKmerIndexer &operator=(ApproxKmerIndexer &&) = delete;

  [[nodiscard]] KmerIndexes extract(const std::vector<Contig> &contigs,
                                    const std::optional<std::vector<Contig>> &readset_optional,
                                    logging::Logger &logger) const {
      const kmer_filter::KmerFilterBuilder kmer_filter_builder
          {nthreads, hasher, common_params, kmer_indexer_params};
      logger.info() << "Creating kmer filter\n";
      const kmer_filter::KmerFilter
          kmer_filter = kmer_filter_builder.GetKmerFilter(contigs, logger);
      logger.info()
          << "Finished creating kmer filter. Using it to build kmer indexes\n";
      KmerIndexes kmer_indexes = GetKmerIndexes(contigs, kmer_filter, logger);
      if (readset_optional.has_value()) {
          // Careful mode
          logger.info()
              << "Careful mode requested. Filtering high multiplicity unique k-mers\n";
          const std::vector<Contig> &readset = readset_optional.value();
          BanHighFreqUniqueKmers(contigs, readset, kmer_indexes, logger);
      }
      return kmer_indexes;
  }
};

}// End namespace veritymap::kmer_index::approx_kmer_indexer