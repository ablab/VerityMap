//
// Created by Andrey Bzikadze on 03/31/21.
//

#pragma once

#include "../rolling_hash.hpp"
#include "kmer_index.hpp"
#include "../config/config.hpp"

namespace tandem_mapper::kmer_index::sketch_contigs {

    template <typename htype>
    class SketchContig {
        const Contig & contig;
        const RollingHash<htype> & hasher;
        size_t k;
        sketch::cm::ccm_t cms;

    public:
        SketchContig(const Contig & contig_,
                     BloomFilter & once_filter,
                     BloomFilter & ban_filter,
                     const RollingHash<htype> & hasher_,
                     const size_t k_,
                     const size_t max_cnt,
                     const double exp_base,
                     const int nhash) : contig{contig_}, hasher{hasher_}, k{k_},
                                        cms{static_cast<int>(ceil(log2(static_cast<double>(max_cnt)))),
                                            static_cast<int>(ceil(log2(std::exp(exp_base) *
                                                                          static_cast<double>(contig_.size())))),
                                            nhash} {
            if (contig.size() < hasher.k) {
                return;
            }
            KWH<htype> kwh(hasher, contig.seq, 0);
            while(true) {
                const htype fhash { kwh.get_fhash() };
                const htype rhash { kwh.get_rhash() };

                ban_filter.insert(rhash);
                if (ban_filter.contains(fhash)) {
                    // skip
                } else if (once_filter.contains(fhash) and cms.est_count(fhash) == 0) {
                    // kmer appeared in some other contig but not in the current one
                    ban_filter.insert(fhash);
                } else {
                    cms.add(fhash);
                    once_filter.insert(fhash);
                }

                if (!kwh.hasNext()) {
                    break;
                }
                kwh = kwh.next();
            }
        }

        [[nodiscard]] const sketch::ccm_t & get_cms() const { return cms; }
    };

    template <typename htype>
    class SketchContigs {
        const std::vector<Contig> & contigs;
        std::vector<SketchContig<htype>> sketch_contigs;
        const RollingHash<htype> & hasher;
        size_t k {0};
        size_t max_cnt {1};
        BloomFilter once_filter;
        BloomFilter ban_filter;
        static BloomParameters _get_filter_params(const std::vector<Contig> & contigs,
                                                  const size_t k, const double fpp) {
            size_t N {0};
            for (const Contig & contig : contigs) {
                N += contig.size() - k + 1;
            }
            BloomParameters params;
            params.projected_element_count = N;
            params.false_positive_probability = fpp;
            params.compute_optimal_parameters();
            return params;
        }

    public:
        SketchContigs(const std::vector<Contig> & contigs_,
                      const RollingHash<htype> & hasher_,
                      const size_t k_,
                      const size_t max_cnt_,
                      const double fpp,
                      const double exp_base,
                      const int nhash): contigs{contigs_}, hasher{hasher_}, k{k_}, max_cnt{max_cnt_},
                                        once_filter{_get_filter_params(contigs_, k_, fpp)},
                                        ban_filter {_get_filter_params(contigs_, k_, fpp)} {
            for (const Contig & contig : contigs) {
                sketch_contigs.emplace_back(contig, once_filter, ban_filter, hasher, k, max_cnt, exp_base, nhash);
            }
        }

        [[nodiscard]] std::vector<tandem_mapper::kmer_index::KmerIndex> get_kmer_indexes(const size_t step_size,
                                                                                         const size_t window_size,
                                                                                         const double uniq_frac) const {
            enum class KmerType { unique, rare, frequent, not_solid };
            struct KmerWindow {
                const size_t length {1};
                size_t tot_uniq {0};
                std::deque<KmerType> deque;

                explicit KmerWindow(const size_t length_): length{length_}, deque{length_, KmerType::not_solid} {
                    VERIFY(length >= 1);
                }

                [[nodiscard]] double get_uniq_frac() const {
                    return static_cast<double>(tot_uniq) / length;
                }

                void popnpush(KmerType kmer_type) {
                    KmerType kmer_type_front { deque.front() };
                    deque.pop_front();
                    if (kmer_type_front == KmerType::unique) {
                        --tot_uniq;
                    }
                    deque.push_back(std::move(kmer_type));
                    if (kmer_type == KmerType::unique) {
                        ++tot_uniq;
                    }
                }
            };

            std::vector<tandem_mapper::kmer_index::KmerIndex> kmer_indexes;
            VERIFY(contigs.size() == sketch_contigs.size());
            for (auto [itcontig, itsc] = std::pair{contigs.cbegin(), sketch_contigs.cbegin()};
                     itcontig != contigs.cend();
                     ++itcontig, ++itsc) {
                tandem_mapper::kmer_index::KmerIndex & kmer_index { kmer_indexes.emplace_back() };
                const Contig & contig = *itcontig;
                const SketchContig<htype> & sketch_contig = *itsc;
                if (contig.size() < hasher.k) {
                    continue;
                }
                KmerWindow kmer_window {window_size};
                KWH<htype> kwh(hasher, contig.seq, 0);
                while(true) {
                    const htype fhash = kwh.get_fhash();
                    const htype rhash = kwh.get_rhash();
                    const KmerType kmer_type = [&fhash, &rhash, &sketch_contig, this]() {
                        if (ban_filter.contains(fhash)) {
                            return KmerType::not_solid;
                        }
                        const size_t fcnt{sketch_contig.get_cms().est_count(fhash)};
                        if (fcnt > max_cnt) {
                            return KmerType::frequent;
                        }
                        return fcnt == 1 ? KmerType::unique : KmerType::rare;
                    }();
                    kmer_window.popnpush(kmer_type);
                    const bool is_solid = (kmer_type == KmerType::unique) or (kmer_type == KmerType::rare);

                    if (is_solid and ((kmer_window.get_uniq_frac() < uniq_frac) or (kwh.pos % step_size == 0))) {
                        kmer_index[fhash].emplace_back(kwh.pos);
                    }

                    if (!kwh.hasNext()) {
                        break;
                    }
                    kwh = kwh.next();
                }
            }
            return kmer_indexes;
        }
    };

    template <typename htype>
    std::vector<KmerIndex> get_rare_kmers_approx(const std::vector<Contig> & contigs,
                                                        const RollingHash<htype> & hasher,
                                                        const Config::CommonParams & common_params,
                                                        const Config::KmerIndexerParams & kmer_indexer_params) {
        const Config::KmerIndexerParams::ApproximateKmerIndexerParams & approximate_kmer_indexer_params =
                kmer_indexer_params.approximate_kmer_indexer_params;
        SketchContigs<htype> sketch_contigs{contigs, hasher, common_params.k, kmer_indexer_params.max_rare_cnt_target,
                                            approximate_kmer_indexer_params.false_positive_probability,
                                            approximate_kmer_indexer_params.exp_base,
                                            approximate_kmer_indexer_params.nhash};
        return sketch_contigs.get_kmer_indexes(kmer_indexer_params.k_step_size,
                                               kmer_indexer_params.k_window_size,
                                               kmer_indexer_params.window_unique_density);
    }

} // End namespace tandem_mapper::kmer_index::sketch_contigs