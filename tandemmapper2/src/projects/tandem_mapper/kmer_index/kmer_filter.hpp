//
// Created by Andrey Bzikadze on 10/26/21.
//

#pragma once

#include "../cms_utils.hpp"

namespace tandem_mapper::kmer_index::kmer_filter {

enum class KmerType { unique,
                      rare,
                      frequent,
                      banned };

class KmerFilter {
  std::vector<BloomFilter> once_filters;
  std::vector<BloomFilter> ban_filters;
  std::vector<std::vector<sketch::cm::ccm_t>> cmss;

  template<typename htype>
  friend class KmerFilterBuilder;

 public:
  KmerFilter(std::vector<BloomFilter> once_filters,
             std::vector<BloomFilter> ban_filters) : once_filters{std::move(once_filters)},
                                                     ban_filters{std::move(ban_filters)} {
    VERIFY(once_filters.size() == ban_filters.size());
    if (not cmss.empty()) {
      VERIFY(ban_filters.size() == cmss.front().size());
    }
  }

  template<typename htype>
  KmerType GetKmerType(const size_t ctg_ind, const htype fhash, const size_t i, const size_t max_rare_cnt) const {
    const BloomFilter &once_filter = once_filters[i];
    const BloomFilter &ban_filter = ban_filters[i];
    const sketch::cm::ccm_t &cms = cmss[ctg_ind][i];
    if (ban_filter.contains(fhash)) {
      return KmerType::banned;
    }
    const size_t fcnt{cms.est_count(fhash)};
    if (fcnt > max_rare_cnt) {
      return KmerType::frequent;
    }
    return fcnt == 1 ? KmerType::unique : KmerType::rare;
  }
};

template<typename htype>
class KmerFilterBuilder {
  size_t nthreads{0};
  const RollingHash<htype> &hasher;
  Config::CommonParams common_params;
  Config::KmerIndexerParams kmer_indexer_params;

 private:
  [[nodiscard]] BloomParameters GetBloomParams(const uint64_t tot_len) const {
    const Config::KmerIndexerParams::ApproximateKmerIndexerParams &approx_kmer_indexer_params =
        kmer_indexer_params.approximate_kmer_indexer_params;
    BloomParameters params;
    params.projected_element_count = tot_len / nthreads;
    params.false_positive_probability = approx_kmer_indexer_params.false_positive_probability;
    params.compute_optimal_parameters();
    return params;
  }

  [[nodiscard]] KmerFilter InitKmerFilter(const std::vector<Contig> &contigs) const {
    uint64_t tot_len{0};
    for (const Contig &contig : contigs) {
      tot_len += contig.size();
    }
    const BloomParameters kBloomParameters = GetBloomParams(tot_len);
    std::vector<BloomFilter> once_filters;
    (nthreads, BloomFilter(kBloomParameters));
    std::vector<BloomFilter> ban_filters;
    (nthreads, BloomFilter(kBloomParameters));
    for (size_t i = 0; i < nthreads; ++i) {
      once_filters.emplace_back(kBloomParameters);
      ban_filters.emplace_back(kBloomParameters);
    }

    return {std::move(once_filters), std::move(ban_filters)};
  }

  void AddContigToFilter(KmerFilter &kmer_filter,
                         const Contig &contig) const {
    const cms_utils::CMSParams kCmsParams(common_params, kmer_indexer_params, contig.size(), nthreads);
    std::vector<sketch::cm::ccm_t> cms;
    for (size_t i = 0; i < nthreads; ++i) {
      cms.emplace_back(kCmsParams.nbits, kCmsParams.l2sz, kCmsParams.nhash);
    }
    kmer_filter.cmss.emplace_back(std::move(cms));

    if (contig.size() < common_params.k) {
      return;
    }

    std::vector<std::vector<std::pair<htype, htype>>> hashes(nthreads);
    std::vector<size_t> sizes(nthreads, 0);

    auto process_chunk = [&kmer_filter, &sizes, &hashes](const size_t i) {
      BloomFilter &ban_f = kmer_filter.ban_filters[i];
      BloomFilter &once_f = kmer_filter.once_filters[i];
      sketch::cm::ccm_t &sketch = kmer_filter.cmss.back()[i];
      const std::vector<std::pair<htype, htype>> &hashes_th = hashes[i];
      const size_t size = sizes[i];
      for (int j = 0; j < size; ++j) {
        const auto [fhash, rhash] = hashes_th[j];
        ban_f.insert(rhash);
        if (ban_f.contains(fhash)) {
          // skip
        } else if (once_f.contains(fhash) and sketch.est_count(fhash) == 0) {
          // kmer appeared in some other contig but not in the current one
          ban_f.insert(fhash);
        } else {
          sketch.add(fhash);
          once_f.insert(fhash);
        }
      }
    };

    const size_t chunk_size = kmer_indexer_params.approximate_kmer_indexer_params.chunk_size;

    KWH<htype> kwh({hasher, contig.seq, 0});
    while (true) {
      for (size_t cnt = 0; cnt < chunk_size; ++cnt) {
        const htype fhash = kwh.get_fhash();
        const htype rhash = kwh.get_rhash();
        const size_t ithread = ((fhash * rhash) % (2 * nthreads)) / 2;
        if (hashes[ithread].size() == sizes[ithread]) {
          hashes[ithread].emplace_back(fhash, rhash);
        } else {
          hashes[ithread][sizes[ithread]] = {fhash, rhash};
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

      if (not kwh.hasNext()) {
        break;
      }

      std::fill(sizes.begin(), sizes.end(), 0);
    }
  }

 public:
  KmerFilterBuilder(size_t nthreads,
                    const RollingHash<htype> &hasher,
                    const Config::CommonParams &common_params,
                    const Config::KmerIndexerParams &kmer_indexer_params) : nthreads(nthreads),
                                                                            hasher(hasher),
                                                                            common_params(common_params),
                                                                            kmer_indexer_params(kmer_indexer_params) {}

  KmerFilterBuilder(const KmerFilterBuilder &) = delete;
  KmerFilterBuilder(KmerFilterBuilder &&) = delete;
  KmerFilterBuilder &operator=(const KmerFilterBuilder &) = delete;
  KmerFilterBuilder &operator=(KmerFilterBuilder &&) = delete;

  [[nodiscard]] KmerFilter GetKmerFilter(const std::vector<Contig> &contigs) const {
    KmerFilter kmer_filter = InitKmerFilter(contigs);
    for (const Contig &contig : contigs) {
      AddContigToFilter(kmer_filter, contig);
    }
    return kmer_filter;
  }
};

}// End namespace tandem_mapper::kmer_index::kmer_filter