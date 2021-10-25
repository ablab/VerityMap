//
// Created by Andrey Bzikadze on 10/19/21.
//

#pragma once

#include "sketch_contigs.hpp"

namespace tandem_mapper::kmer_index::sketch_contigs {
    template<typename htype>
    class SketchContig;
}

namespace tandem_mapper::kmer_index::kmer_type {

    enum class KmerType { unique, rare, frequent, banned };

    template<typename htype>
    KmerType get_kmer_type(
            const htype fhash,
            const sketch::cm::ccm_t & cms,
            const BloomFilter & ban_filter,
            const size_t min_freq_cnt) {
        if (ban_filter.contains(fhash)) {
            return KmerType::banned;
        }
        const size_t fcnt{cms.est_count(fhash)};
        if (fcnt > min_freq_cnt) {
            return KmerType::frequent;
        }
        return fcnt == 1 ? KmerType::unique : KmerType::rare;
    }

    template<typename htype>
    KmerType get_kmer_type(
            const htype fhash,
            const std::vector<tandem_mapper::kmer_index::sketch_contigs::SketchContig<htype>> & sketch_contigs,
            const BloomFilter & ban_filter,
            const size_t min_freq_cnt) {
        if (ban_filter.contains(fhash)) {
            return KmerType::banned;
        }

        size_t fcnt {0};
        for (const tandem_mapper::kmer_index::sketch_contigs::SketchContig<htype> & sketch_contig : sketch_contigs) {
            fcnt = std::max(fcnt, sketch_contig.get_cms().est_count(fhash));
        }
        if (fcnt > min_freq_cnt) {
            return KmerType::frequent;
        }
        return fcnt == 1 ? KmerType::unique : KmerType::rare;
    }


} // End namespace tandem_mapper::kmer_index::kmer_type