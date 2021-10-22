//
// Created by Andrey Bzikadze on 03/31/21.
//

#pragma once
#include <omp.h>
#include <ctime>

#include "kmer_index.hpp"
#include "kmer_type.hpp"
#include "kmer_window.hpp"
#include "../rolling_hash.hpp"
#include "../config/config.hpp"
#include "common/coverage_utils.hpp"

namespace tandem_mapper::kmer_index::sketch_contigs {

    template <typename htype>
    class SketchContig {
        const Contig & contig;
        const RollingHash<htype> & hasher;
        size_t k;
        std::vector<sketch::cm::ccm_t> cms;

    public:
        SketchContig(const Contig & contig_,
                     std::vector<BloomFilter> & once_filter,
                     std::vector<BloomFilter> & ban_filter,
                     const RollingHash<htype> & hasher_,
                     const size_t k_,
                     const size_t max_cnt,
                     const double exp_base,
                     const int nhash,
                     const int nthreads,
                     const size_t chunk_size) :
                         contig{contig_}, hasher{hasher_}, k{k_},
                         cms(nthreads,
                             sketch::cm::ccm_t(static_cast<int>(ceil(log2(static_cast<double>(max_cnt)))),
                                               static_cast<int>(ceil(log2(std::exp(exp_base) *
                                                           static_cast<double>(contig_.size()) / nthreads))),
                                                           nhash)
                             ) {
            if (contig.size() < hasher.k) {
                return;
            }

            std::vector<std::vector<std::pair<htype, htype>>> hashes(nthreads);
            std::vector<size_t> sizes(nthreads, 0);

            auto process_chunk = [&ban_filter, &once_filter, &sizes, &hashes, this](size_t i) {
                BloomFilter & ban_f = ban_filter[i];
                BloomFilter & once_f = once_filter[i];
                sketch::cm::ccm_t & sketch = cms[i];
                std::vector<std::pair<htype, htype>> & hashes_th = hashes[i];
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

            const size_t thread_chunk_size = std::max<size_t>(1, chunk_size / nthreads);

            KWH<htype> kwh({hasher, contig.seq, 0});
            while (true) {
                for (size_t cnt = 0; cnt < chunk_size; ++cnt) {
                    const htype fhash = kwh.get_fhash();
                    const htype rhash = kwh.get_rhash();
                    const size_t ithread = fhash % nthreads;
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
                for (auto & thread : threads) {
                    thread.join();
                }

                if (not kwh.hasNext()) {
                    break;
                }

                std::fill(sizes.begin(), sizes.end(), 0);
            }
        }

        [[nodiscard]] const std::vector<sketch::cm::ccm_t> & get_cms() const {
            return cms;
        }
    };


    template <typename htype>
    class SketchContigs {
        const std::vector<Contig> & contigs;
        std::vector<SketchContig<htype>> sketch_contigs;
        const RollingHash<htype> & hasher;
        size_t k {0};
        size_t max_cnt {1};
        std::vector<BloomFilter> once_filter;
        std::vector<BloomFilter> ban_filter;
        double careful_upper_bnd_cov_mult;

        static BloomParameters _get_filter_params(const std::vector<Contig> & contigs,
                                                  const size_t k, const double fpp, const uint64_t nthreads) {
            size_t N {0};
            for (const Contig & contig : contigs) {
                N += contig.size() - k + 1;
            }
            BloomParameters params;
            params.projected_element_count = N / nthreads;
            params.false_positive_probability = fpp;
            params.compute_optimal_parameters();
            return params;
        }

    public:
        SketchContigs(const std::vector<Contig> & contigs,
                      const std::vector<Contig> & readset,
                      const RollingHash<htype> & hasher,
                      const size_t k,
                      const size_t max_cnt,
                      const double fpp,
                      const double exp_base,
                      const int nhash,
                      const double careful_upper_bnd_cov_mult,
                      const uint64_t nthreads,
                      const size_t chunk_size): contigs{contigs}, hasher{hasher}, k{k}, max_cnt{max_cnt},
                                                once_filter(nthreads, _get_filter_params(contigs, k, fpp, nthreads)),
                                                ban_filter (nthreads, _get_filter_params(contigs, k, fpp, nthreads)),
                                                careful_upper_bnd_cov_mult {careful_upper_bnd_cov_mult} {
            for (const Contig & contig : contigs) {
                sketch_contigs.emplace_back(contig, once_filter, ban_filter, hasher, k, max_cnt, exp_base, nhash, nthreads, chunk_size);
            }
            // ban_high_freq_unique_kmers(contigs, readset, exp_base, nhash, nthreads);
        }

        uint64_t get_n_unique_kmers() {
            using namespace tandem_mapper::kmer_index::kmer_type;
            uint64_t n_unique_kmers { 0 };
            for (auto [itcontig, itsc] = std::pair{contigs.cbegin(), sketch_contigs.cbegin()};
                 itcontig != contigs.cend();
                 ++itcontig, ++itsc) {
                const Contig & contig = *itcontig;
                const SketchContig<htype> & sketch_contig = *itsc;
                if (contig.size() < hasher.k) {
                    continue;
                }
                KWH<htype> kwh(hasher, contig.seq, 0);
                while(true) {
                    const htype fhash = kwh.get_fhash();
                    const htype rhash = kwh.get_rhash();

                    const KmerType kmer_type = get_kmer_type(fhash, rhash, sketch_contig, ban_filter, max_cnt);
                    if (kmer_type == KmerType::unique) {
                        ++n_unique_kmers;
                    }

                    if (!kwh.hasNext()) {
                        break;
                    }
                    kwh = kwh.next();
                }
            }
            return n_unique_kmers;
        }

        [[nodiscard]] tandem_mapper::kmer_index::KmerIndexes get_kmer_indexes(const std::vector<Contig> & readset,
                                                                              const size_t nthreads,
                                                                              const size_t step_size,
                                                                              const size_t window_size,
                                                                              const double uniq_frac,
                                                                              const size_t chunk_size) const {
            using tandem_mapper::kmer_index::kmer_window::KmerWindow;
            using tandem_mapper::kmer_index::KmerIndex;
            tandem_mapper::kmer_index::KmerIndexes kmer_indexes_all;
            VERIFY(contigs.size() == sketch_contigs.size());
            for (auto [itcontig, itsc] = std::pair{contigs.cbegin(), sketch_contigs.cbegin()};
                     itcontig != contigs.cend();
                     ++itcontig, ++itsc) {
                std::vector<KmerIndex> kmer_indexes(nthreads);
                const Contig & contig = *itcontig;
                const SketchContig<htype> & sketch_contig = *itsc;
                const std::vector<sketch::cm::ccm_t> & cms = sketch_contig.get_cms();
                if (contig.size() < hasher.k) {
                    continue;
                }

                std::vector<std::vector<std::tuple<htype, htype, size_t>>> hashes_pos(nthreads);
                std::vector<size_t> sizes(nthreads, 0);
                std::vector<KmerWindow> kmer_windows(nthreads, KmerWindow(window_size, nthreads));

                auto process_chunk = [this,
                                      &cms,
                                      &sizes,
                                      &hashes_pos,
                                      &kmer_windows,
                                      &kmer_indexes,
                                      &nthreads,
                                      &uniq_frac,
                                      &step_size](size_t i) {
                    const BloomFilter & ban_f = ban_filter[i];
                    const BloomFilter & once_f = once_filter[i];
                    const sketch::cm::ccm_t & sketch = cms[i];
                    const std::vector<std::tuple<htype, htype, size_t>> & hashes_pos_th = hashes_pos[i];
                    KmerWindow & kmer_window = kmer_windows[i];
                    KmerIndex & kmer_index = kmer_indexes[i];
                    const size_t size = sizes[i];
                    for (int j = 0; j < size; ++j) {
                        const auto[fhash, rhash, pos] = hashes_pos_th[j];
                        const size_t ithread = fhash % nthreads;
                        const kmer_type::KmerType kmer_type =
                                tandem_mapper::kmer_index::kmer_type::get_kmer_type(fhash,
                                                                                    rhash,
                                                                                    cms[ithread],
                                                                                    ban_filter[ithread],
                                                                                    max_cnt);
                        kmer_window.popnpush(kmer_type);
                        const bool is_solid = (kmer_type == kmer_type::KmerType::unique) or
                                              (kmer_type == kmer_type::KmerType::rare);

                        if (is_solid and ((kmer_window.get_uniq_frac() < uniq_frac) or (pos % step_size == 0))) {
                            kmer_index[fhash].emplace_back(pos);
                        }
                    }
                };

                const size_t thread_chunk_size = std::max<size_t>(1, chunk_size / nthreads);

                KWH<htype> kwh({hasher, contig.seq, 0});
                while (true) {
                    for (size_t cnt = 0; cnt < chunk_size; ++cnt) {
                        const htype fhash = kwh.get_fhash();
                        const htype rhash = kwh.get_rhash();
                        const size_t ithread = fhash % nthreads;
                        if (hashes_pos[ithread].size() == sizes[ithread]) {
                            hashes_pos[ithread].emplace_back(fhash, rhash, kwh.pos);
                        } else {
                            hashes_pos[ithread][sizes[ithread]] = {fhash, rhash, kwh.pos};
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
                    for (auto & thread : threads) {
                        thread.join();
                    }

                    if (not kwh.hasNext()) {
                        break;
                    }

                    std::fill(sizes.begin(), sizes.end(), 0);
                }

                KmerIndex kmer_index_merged;
                for (auto & kmer_index : kmer_indexes) {
                    kmer_index_merged.merge(kmer_index);
                }
                kmer_indexes_all.emplace_back(std::move(kmer_index_merged));
            }
            return kmer_indexes_all;
        }
    };

    template <typename htype>
    std::vector<KmerIndex> get_rare_kmers_approx(const std::vector<Contig> & contigs,
                                                 const std::vector<Contig> & readset,
                                                 const size_t nthreads,
                                                 const RollingHash<htype> & hasher,
                                                 const Config::CommonParams & common_params,
                                                 const Config::KmerIndexerParams & kmer_indexer_params) {
        const Config::KmerIndexerParams::ApproximateKmerIndexerParams & approximate_kmer_indexer_params =
                kmer_indexer_params.approximate_kmer_indexer_params;
        SketchContigs<htype> sketch_contigs{contigs, readset, hasher,
                                            common_params.k,
                                            kmer_indexer_params.max_rare_cnt_target,
                                            approximate_kmer_indexer_params.false_positive_probability,
                                            approximate_kmer_indexer_params.exp_base,
                                            approximate_kmer_indexer_params.nhash,
                                            kmer_indexer_params.careful_upper_bnd_cov_mult,
                                            nthreads,
                                            approximate_kmer_indexer_params.chunk_size};
        std::vector<tandem_mapper::kmer_index::KmerIndex> kmer_indexes =
                sketch_contigs.get_kmer_indexes(readset, nthreads, kmer_indexer_params.k_step_size,
                                                kmer_indexer_params.k_window_size,
                                                kmer_indexer_params.window_unique_density,
                                                approximate_kmer_indexer_params.chunk_size);

        return kmer_indexes;
    }

} // End namespace tandem_mapper::kmer_index::sketch_contigs





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
