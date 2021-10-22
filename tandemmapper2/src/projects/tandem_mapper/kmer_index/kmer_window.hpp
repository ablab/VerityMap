//
// Created by Andrey Bzikadze on 10/19/21.
//


#pragma once

namespace tandem_mapper::kmer_index::kmer_window {
    struct KmerWindow {
        size_t tot_uniq {0};
        std::deque<kmer_type::KmerType> deque;

        explicit KmerWindow(const size_t length, const size_t nthreads):
                deque{length / nthreads, kmer_type::KmerType::banned} {
            VERIFY(length >= 1);
            VERIFY(nthreads >= 1);
        }

        [[nodiscard]] double get_uniq_frac() const {
            return static_cast<double>(tot_uniq) / deque.size();
        }

        void popnpush(kmer_type::KmerType kmer_type) {
            kmer_type::KmerType kmer_type_front { deque.front() };
            deque.pop_front();
            if (kmer_type_front == kmer_type::KmerType::unique) {
                --tot_uniq;
            }
            deque.push_back(kmer_type);
            if (kmer_type == kmer_type::KmerType::unique) {
                ++tot_uniq;
            }
        }
    };
} // End namespace tandem_mapper::kmer_index::kmer_window