//
// Created by Andrey Bzikadze on 10/19/21.
//

#pragma once

namespace veritymap::kmer_index::kmer_window {
template<typename htype>
class KmerWindow {
  int64_t half_length{0};
  int64_t tot_regular{0};
  const std::vector<std::tuple<size_t, htype, bool>> &pos_hash_regular;
  typename std::vector<std::tuple<size_t, htype, bool>>::const_iterator next_it;
  typename std::vector<std::tuple<size_t, htype, bool>>::const_iterator cur_it;
  std::deque<std::pair<int64_t, bool>> deque;

  void IncRight() {
    for (const size_t pos = std::get<0>(*cur_it); next_it != pos_hash_regular.cend(); ++next_it) {
      auto [next_pos, _, is_next_regular] = *next_it;
      if (next_pos - pos > half_length) {
        break;
      }
      deque.emplace_back(next_pos, is_next_regular);
      tot_regular += is_next_regular;
      ++next_pos;
    }
  }

  void IncLeft() {
    const size_t pos = std::get<0>(*cur_it);
    while (not deque.empty()) {
      const auto [pos_front, is_regular_front]{deque.front()};
      if (pos - pos_front <= half_length) {
        break;
      }
      deque.pop_front();
      tot_regular -= is_regular_front;
    }
  }

 public:
  KmerWindow(const size_t length, const std::vector<std::tuple<size_t, htype, bool>> &pos_hash_regular)
      : half_length{(int64_t) length / 2},
        pos_hash_regular{pos_hash_regular},
        next_it{pos_hash_regular.cbegin()},
        cur_it{pos_hash_regular.cbegin()} {
    VERIFY(length >= 1);
  }

  [[nodiscard]] double RegularFrac() const { return tot_regular / double(half_length * 2); }

  void Inc() {
    ++cur_it;
    IncLeft();
    IncRight();
  }

  void Reset() {
    next_it = pos_hash_regular.cbegin();
    cur_it = pos_hash_regular.cbegin();
    IncRight();
  }
};
}// End namespace veritymap::kmer_index::kmer_window