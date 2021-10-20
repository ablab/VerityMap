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
        double careful_upper_bnd_cov_mult;

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

        void ban_high_freq_unique_kmers(const std::vector<Contig> & contigs_,
                                        const std::vector<Contig> & readset,
                                        const double exp_base,
                                        const int nhash,
                                        const uint32_t nthreads) {
            // If read-set is not empty, we additionally ban unique k-mers in assembly that have unusually high coverage
            if (readset.empty())
                return;

            uint64_t n_unique_kmers { get_n_unique_kmers() };

            const double coverage { tools::common::coverage_utils::get_coverage(contigs_, readset) };

            const uint max_read_freq = std::max(1., ceil(careful_upper_bnd_cov_mult * coverage));
            const int nbits = std::max(1., ceil(log2(max_read_freq)));
            const int l2sz = ceil(log2(
                    std::exp(exp_base) * ((double) n_unique_kmers)
            ));

            sketch::cm::ccm_t cms {nbits, l2sz, nhash};

            for (const Contig & contig : readset) {
                if (contig.size() < hasher.k) {
                    continue;
                }
                KWH<htype> kwh(hasher, contig.seq, 0);
                while(true) {
                    const htype fhash = kwh.get_fhash();
                    const htype rhash = kwh.get_rhash();
                    std::vector<std::pair<htype, htype>> hashes { { fhash, rhash }, { rhash, fhash } };
                    for (const auto [x, y] : hashes) {
                        kmer_type::KmerType kmer_type =
                                tandem_mapper::kmer_index::kmer_type::get_kmer_type(x, y,
                                                                                    sketch_contigs,
                                                                                    ban_filter,
                                                                                    max_cnt);
                        if (kmer_type == kmer_type::KmerType::unique) {
                            if (ban_filter.contains((x))) {
                                continue;
                            } else {
                                cms.add(x);
                                if (cms.est_count(x) == max_read_freq) {
                                    ban_filter.insert(x);
                                }
                            }
                        }
                    }
                    if (!kwh.hasNext()) {
                        break;
                    }
                    kwh = kwh.next();
                }
            }
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
                      const uint64_t nthreads): contigs{contigs}, hasher{hasher}, k{k}, max_cnt{max_cnt},
                                                once_filter{_get_filter_params(contigs, k, fpp)},
                                                ban_filter {_get_filter_params(contigs, k, fpp)},
                                                careful_upper_bnd_cov_mult {careful_upper_bnd_cov_mult} {
            for (const Contig & contig : contigs) {
                sketch_contigs.emplace_back(contig, once_filter, ban_filter, hasher, k, max_cnt, exp_base, nhash);
            }
            ban_high_freq_unique_kmers(contigs, readset, exp_base, nhash, nthreads);
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
                                                                              const double uniq_frac) const {
            tandem_mapper::kmer_index::KmerIndexes kmer_indexes;
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
                tandem_mapper::kmer_index::kmer_window::KmerWindow kmer_window {window_size};
                KWH<htype> kwh(hasher, contig.seq, 0);
                while(true) {
                    const htype fhash = kwh.get_fhash();
                    const htype rhash = kwh.get_rhash();

                    const kmer_type::KmerType kmer_type =
                        tandem_mapper::kmer_index::kmer_type::get_kmer_type(fhash,
                                                                            rhash,
                                                                            sketch_contig,
                                                                            ban_filter,
                                                                            max_cnt);
                    kmer_window.popnpush(kmer_type);
                    const bool is_solid = (kmer_type == kmer_type::KmerType::unique) or
                            (kmer_type == kmer_type::KmerType::rare);

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
                                            nthreads};
        std::vector<tandem_mapper::kmer_index::KmerIndex> kmer_indexes =
                sketch_contigs.get_kmer_indexes(readset, nthreads, kmer_indexer_params.k_step_size,
                                                kmer_indexer_params.k_window_size,
                                                kmer_indexer_params.window_unique_density);

        return kmer_indexes;
    }

} // End namespace tandem_mapper::kmer_index::sketch_contigs