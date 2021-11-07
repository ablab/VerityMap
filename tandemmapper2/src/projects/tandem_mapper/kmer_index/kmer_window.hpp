//
// Created by Andrey Bzikadze on 10/19/21.
//

#pragma once

namespace tandem_mapper::kmer_index::kmer_window {
struct KmerWindow {
  std::size_t length{0};
  int64_t tot_uniq;
  std::deque<std::pair<size_t, bool>> deque;

  explicit KmerWindow(const size_t length) : length{length},
                                             tot_uniq{0} {
    VERIFY(length >= 1);
  }

  [[nodiscard]] double unique_frac() const {
    return tot_uniq / double(length);
  }

  void add(size_t pos, bool is_unique) {
    while (not deque.empty()) {
      const auto [pos_front, is_unique_front]{deque.front()};
      if (pos - pos_front < length) {
        break;
      }
      deque.pop_front();
      if (is_unique_front) {
        --tot_uniq;
      }
    }
    deque.emplace_back(pos, is_unique);
    if (is_unique) {
      ++tot_uniq;
    }
  }
};
}// End namespace tandem_mapper::kmer_index::kmer_window