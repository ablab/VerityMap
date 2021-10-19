//
// Created by Andrey Bzikadze on 10/18/21.
//

#pragma once

#include "bloom/bloom.hpp"

namespace tandem_mapper::kmer_index::filter_rep_kmers {

    template<typename htype>
    BloomFilter get_bloom_rep_kmers(const Sequence & sequence,
                                    const RollingHash<htype> & hasher,
                                    const double false_positive_probability) {
        if (sequence.size() < hasher.k) {
            return {};
        }
        BloomParameters bloom_params;
        bloom_params.projected_element_count = sequence.size();
        bloom_params.false_positive_probability = false_positive_probability;
        bloom_params.compute_optimal_parameters();

        BloomFilter once_filter{bloom_params};
        BloomFilter twice_filter{bloom_params};

        KWH<htype> kwh(hasher, sequence, 0);
        while (true) {
            const htype hash = kwh.get_fhash();
            if (once_filter.contains(hash)) {
                twice_filter.insert(hash);
            } else {
                once_filter.insert(hash);
            }
            if (not kwh.hasNext()) {
                break;
            }
            kwh = kwh.next();
        }
        return twice_filter;
    }

} // End namespace tandem_mapper::kmer_index::filter_rep_kmers