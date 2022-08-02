//
// Created by Andrey Bzikadze on 06/02/22.
//

#pragma once

#include "../cms_utils.hpp"
#include "sketch/ccm.h"

namespace veritymap::kmer_index::kmer_filter_canon {

enum class KmerTypeCanon { unique, duplicate, rare, frequent };

bool IsKmerTypeCanonRegular(const KmerTypeCanon kmer_type, const bool diploid = false) {
  return kmer_type == KmerTypeCanon::unique or (diploid and kmer_type == KmerTypeCanon::duplicate);
}

class KmerFilterCanon {
  std::vector<std::vector<sketch::cm::ccm_t>> cmss;

  template<typename htype>
  friend class KmerFilterCanonBuilder;

 public:
  template<typename htype>
  [[nodiscard]] KmerTypeCanon GetKmerType(const size_t ctg_ind, const htype canhash, const size_t i,
                                          const size_t max_rare_cnt) const {
    const sketch::cm::ccm_t &cms = cmss[ctg_ind][i];
    const size_t fcnt{cms.est_count(canhash)};
    if (fcnt > max_rare_cnt)
      return KmerTypeCanon::frequent;
    if (fcnt == 1)
      return KmerTypeCanon::unique;
    if (fcnt == 2)
      return KmerTypeCanon::duplicate;
    return KmerTypeCanon::rare;
  }
};

template<typename htype>
class KmerFilterCanonBuilder {
  size_t nthreads{0};
  const RollingHash<htype> &hasher;
  Config::CommonParams common_params;
  Config::KmerIndexerParams kmer_indexer_params;

  void AddContigToFilter(KmerFilterCanon &kmer_filter, const Contig &contig, logging::Logger &logger) const {
    const cms_utils::CMSParams kCmsParams(common_params, kmer_indexer_params, contig.size(), nthreads);
    std::vector<sketch::cm::ccm_t> cms;
    for (size_t i = 0; i < nthreads; ++i) { cms.emplace_back(kCmsParams.nbits, kCmsParams.l2sz, kCmsParams.nhash); }
    kmer_filter.cmss.emplace_back(std::move(cms));

    if (contig.size() < common_params.k) {
      return;
    }

    std::vector<std::vector<htype>> hashes(nthreads);
    std::vector<size_t> sizes(nthreads, 0);

    auto process_chunk = [&kmer_filter, &sizes, &hashes](const size_t i) {
      sketch::cm::ccm_t &sketch = kmer_filter.cmss.back()[i];
      const std::vector<htype> &hashes_th = hashes[i];
      const size_t size = sizes[i];
      for (int j = 0; j < size; ++j) {
        const htype hash = hashes_th[j];
        sketch.add(hash);
      }
    };

    const size_t chunk_size = kmer_indexer_params.approximate_kmer_indexer_params.chunk_size;

    KWH<htype> kwh({hasher, contig.seq, 0});
    while (true) {
      logger.info() << "Generating task list for chunk starting at pos " << kwh.pos << "\n";
      for (size_t cnt = 0; cnt < chunk_size; ++cnt) {
        const htype hash = kwh.hash();
        const size_t ithread = hash % nthreads;
        if (hashes[ithread].size() == sizes[ithread]) {
          hashes[ithread].emplace_back(hash);
        } else {
          hashes[ithread][sizes[ithread]] = hash;
        }
        ++sizes[ithread];
        if (not kwh.hasNext()) {
          break;
        }
        kwh = kwh.next();
      }
      logger.info() << "Parallel run for chunk\n";
      std::vector<std::thread> threads(nthreads);
      for (size_t i = 0; i < threads.size(); ++i) { threads[i] = std::thread(process_chunk, i); }
      for (auto &thread : threads) { thread.join(); }

      if (not kwh.hasNext()) {
        break;
      }

      std::fill(sizes.begin(), sizes.end(), 0);
    }
  }

 public:
  KmerFilterCanonBuilder(size_t nthreads, const RollingHash<htype> &hasher, const Config::CommonParams &common_params,
                         const Config::KmerIndexerParams &kmer_indexer_params)
      : nthreads(nthreads),
        hasher(hasher),
        common_params(common_params),
        kmer_indexer_params(kmer_indexer_params) {}

  KmerFilterCanonBuilder(const KmerFilterCanonBuilder &) = delete;
  KmerFilterCanonBuilder(KmerFilterCanonBuilder &&) = delete;
  KmerFilterCanonBuilder &operator=(const KmerFilterCanonBuilder &) = delete;
  KmerFilterCanonBuilder &operator=(KmerFilterCanonBuilder &&) = delete;

  [[nodiscard]] KmerFilterCanon GetKmerFilterCanon(const std::vector<Contig> &contigs, logging::Logger &logger) const {
    logger.info() << "Init filter\n";
    KmerFilterCanon kmer_filter;
    logger.info() << "Start adding contigs to filter\n";
    for (const Contig &contig : contigs) {
      logger.info() << "Add contig " << contig.id << "\n";
      AddContigToFilter(kmer_filter, contig, logger);
    }
    return kmer_filter;
  }
};

}// namespace veritymap::kmer_index::kmer_filter_canon