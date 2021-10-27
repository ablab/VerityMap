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

namespace tandem_mapper::kmer_index::sketch_contigs {

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
                                       KWH<htype> &kwh) const {
    std::vector<size_t> sizes(nthreads, 0);
    auto process_chunk = [this, &hashes_pos, &sizes, &kmer_filter](size_t i) {
      std::vector<HashPosType> &hashes_pos_th = hashes_pos[i];
      const size_t size = sizes[i];
      for (int j = 0; j < size; ++j) {
        const htype fhash = hashes_pos_th[j].fhash;
        hashes_pos_th[j].kmer_type = kmer_filter.GetKmerType(fhash, i, kmer_indexer_params.max_rare_cnt_target);
      }
    };

    const size_t chunk_size = kmer_indexer_params.approximate_kmer_indexer_params.chunk_size;

    for (size_t cnt = 0; cnt < chunk_size; ++cnt) {
      const htype fhash = kwh.get_fhash();
      const htype rhash = kwh.get_rhash();
      // TODO currently only half of threads is being used
      const size_t ithread = (fhash * rhash) % nthreads;
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

    for (int i = 0; i < nthreads; ++i) {
      std::cout << i << " " << sizes[i] << "\n";
    }
    return sizes;
  }

  [[nodiscard]] KmerIndex GetKmerIndex(const Contig &contig,
                                       const kmer_filter::KmerFilter &kmer_filter) const {
    if (contig.size() < hasher.k) {
      return {};
    }

    std::vector<std::vector<HashPosType>> hashes_pos(nthreads);

    KmerIndex kmer_index;
    KWH<htype> kwh({hasher, contig.seq, 0});
    const size_t window_size = kmer_indexer_params.k_window_size;
    kmer_window::KmerWindow kmer_window(window_size);
    while (true) {
      std::vector<size_t> sizes = BinHashesInChunk(hashes_pos, kmer_filter, kwh);
      std::cout << kwh.pos << "\n";

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

      std::sort(pos_hash_uniq.begin(), pos_hash_uniq.end());

      for (const auto &[pos, hash, is_unique] : pos_hash_uniq) {
        kmer_window.add(pos, is_unique);
        if ((kmer_window.unique_frac() < kmer_indexer_params.window_unique_density) or (pos % window_size == 0)) {
          kmer_index[hash].emplace_back(pos);
        }
      }

      if (not kwh.hasNext()) {
        break;
      }
    }
    return kmer_index;
  }

  [[nodiscard]] KmerIndexes GetKmerIndexes(const std::vector<Contig> &contigs,
                                           const kmer_filter::KmerFilter &kmer_filter) const {
    KmerIndexes kmer_indexes;
    for (const Contig &contig : contigs) {
      kmer_indexes.emplace_back(GetKmerIndex(contig, kmer_filter));
    }
    return kmer_indexes;
  }

 public:
  ApproxKmerIndexer(const size_t nthreads,
                    const RollingHash<htype> &hasher,
                    const Config::CommonParams &common_params,
                    const Config::KmerIndexerParams &kmer_indexer_params) : nthreads{nthreads},
                                                                            hasher{hasher},
                                                                            common_params{common_params},
                                                                            kmer_indexer_params{
                                                                                kmer_indexer_params} {}

  ApproxKmerIndexer(const ApproxKmerIndexer &) = delete;
  ApproxKmerIndexer(ApproxKmerIndexer &&) = delete;
  ApproxKmerIndexer &operator=(const ApproxKmerIndexer &) = delete;
  ApproxKmerIndexer &operator=(ApproxKmerIndexer &&) = delete;

  // TODO add careful mode
  // TODO change readset to optional
  [[nodiscard]] KmerIndexes extract(const std::vector<Contig> &contigs,
                                    const std::vector<Contig> &readset) const {
    const kmer_filter::KmerFilterBuilder kmer_filter_builder{nthreads, hasher, common_params, kmer_indexer_params};
    const kmer_filter::KmerFilter kmer_filter = kmer_filter_builder.GetKmerFilter(contigs);
    return GetKmerIndexes(contigs, kmer_filter);
  }
};

}// End namespace tandem_mapper::kmer_index::sketch_contigs

// uint64_t get_n_unique_kmers() {
//   using namespace tandem_mapper::kmer_index::kmer_filter;
//   uint64_t n_unique_kmers{0};
//   for (auto [itcontig, itsc] = std::pair{contigs.cbegin(), sketch_contigs.cbegin()};
//        itcontig != contigs.cend();
//        ++itcontig, ++itsc) {
//     const Contig &contig = *itcontig;
//     const SketchContig<htype> &sketch_contig = *itsc;
//     if (contig.size() < hasher.k) {
//       continue;
//     }
//     KWH<htype> kwh(hasher, contig.seq, 0);
//     while (true) {
//       const htype fhash = kwh.get_fhash();
//       const htype rhash = kwh.get_rhash();
//
//       const KmerType kmer_type = get_kmer_type(fhash, rhash, sketch_contig, ban_filter, max_cnt);
//       if (kmer_type == KmerType::unique) {
//         ++n_unique_kmers;
//       }
//
//       if (!kwh.hasNext()) {
//         break;
//       }
//       kwh = kwh.next();
//     }
//   }
//   return n_unique_kmers;
// }
// void ban_high_freq_unique_kmers(const std::vector<Contig> & contigs_,
//                                 const std::vector<Contig> & readset,
//                                 const double exp_base,
//                                 const int nhash,
//                                 const uint32_t nthreads) {
//     // If read-set is not empty, we additionally ban unique k-mers in assembly that have unusually high coverage
//     if (readset.empty())
//         return;

//     uint64_t n_unique_kmers { get_n_unique_kmers() };

//     const double coverage { tools::common::coverage_utils::get_coverage(contigs_, readset) };

//     const uint max_read_freq = std::max(1., ceil(careful_upper_bnd_cov_mult * coverage));
//     const int nbits = std::max(1., ceil(log2(max_read_freq)));
//     const int l2sz = ceil(log2(
//             std::exp(exp_base) * ((double) n_unique_kmers)
//     ));

//     sketch::cm::ccm_t cms {nbits, l2sz, nhash};

//     for (const Contig & contig : readset) {
//         if (contig.size() < hasher.k) {
//             continue;
//         }
//         KWH<htype> kwh(hasher, contig.seq, 0);
//         while(true) {
//             const htype fhash = kwh.get_fhash();
//             const htype rhash = kwh.get_rhash();
//             std::vector<std::pair<htype, htype>> hashes { { fhash, rhash }, { rhash, fhash } };
//             for (const auto [x, y] : hashes) {
//                 kmer_type::KmerType kmer_type =
//                         tandem_mapper::kmer_index::kmer_type::get_kmer_type(x, y,
//                                                                             sketch_contigs,
//                                                                             ban_filter,
//                                                                             max_cnt);
//                 if (kmer_type == kmer_type::KmerType::unique) {
//                     if (ban_filter.contains((x))) {
//                         continue;
//                     } else {
//                         cms.add(x);
//                         if (cms.est_count(x) == max_read_freq) {
//                             ban_filter.insert(x);
//                         }
//                     }
//                 }
//             }
//             if (!kwh.hasNext()) {
//                 break;
//             }
//             kwh = kwh.next();
//         }
//     }
// }
