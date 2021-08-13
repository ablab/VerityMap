//
// Created by Andrey Bzikadze on 2/22/21.
//

# pragma once

#include "strand.hpp"
#include "config/config.hpp"

namespace tandem_mapper::matches {

    struct Match {
        Config::ChainingParams::match_pos_type target_pos {0};
        int32_t query_pos {0};
        uint8_t target_freq {0};   // TODO Change to "is_unique"

        [[nodiscard]] bool is_unique() const {
            return target_freq == 1;
        }
    };

    inline bool operator< (const Match & lhs, const Match & rhs) {
        return lhs.target_pos < rhs.target_pos;
    }

    inline std::ostream & operator<<(std::ostream & os, const Match & match) {
        os << match.query_pos << "\t" << match.target_pos << "\t" << static_cast<uint16_t>(match.target_freq) << "\n";
        return os;
    }

    using Matches = std::vector<Match>;

    inline std::ostream & operator<<(std::ostream & os, const Matches & matches) {
        size_t prev_pos = 0;
        for (const auto & match : matches) {
            if (match.target_pos - prev_pos > 10) { // TODO: fix to k or k/2
                os << match;
                prev_pos = match.target_pos;
            }
        }
        return os;
    }

    Matches get_matches(const Contig & target,
                        const kmer_index::KmerIndex & target_kmer_index,
                        const Contig & query,
                        const dna_strand::Strand & query_strand,
                        const RollingHash<typename Config::HashParams::htype> & hasher,
                        const size_t max_rare_cnt_query) {
        kmer_index::KmerIndex query_kmer_index_all = query_strand == dna_strand::Strand::forward ?
                kmer_index::get_rare_kmers(query.seq, hasher, max_rare_cnt_query) :
                kmer_index::get_rare_kmers(query.RC().seq, hasher, max_rare_cnt_query);
        Matches matches;
        for (auto && [hash, qpos] : query_kmer_index_all) {
            if (target_kmer_index.contains(hash)) {
                const size_t tf64 = target_kmer_index.at(hash).size();
                VERIFY(tf64 <= std::numeric_limits<uint8_t>::max());
                const auto target_freq = static_cast<uint8_t>(tf64);
                for (const size_t tp : target_kmer_index.at(hash)) {
                    for (const size_t qp : qpos) {
                        VERIFY(qp < std::numeric_limits<int32_t>::max());
                        matches.push_back({ static_cast<Config::ChainingParams::match_pos_type>(tp),
                                            static_cast<int32_t>(qp),
                                            target_freq });
                    }
                }
            }
        }
        std::sort(matches.begin(), matches.end());
        return matches;
    }

} // End namespace tandem_mapper::matches