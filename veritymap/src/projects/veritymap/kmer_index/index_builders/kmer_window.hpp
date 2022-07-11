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

template<class T, class Compare = std::less<T>>
class MinQueue {
  // elements in vector with deque indixes are non-increasing
  std::deque<std::pair<int64_t, T>> deque;
  Compare value_compare;
  int64_t first{0};

 public:
  void PushBack(const int64_t next_pos, const T &next_el) {
    while (not deque.empty() and value_compare(next_el, deque.back().second)) { deque.pop_back(); }
    deque.emplace_back(next_pos, next_el);
  }

  [[nodiscard]] int64_t GetMinIndex() const { return deque.front().first; }
  [[nodiscard]] T GetMin() const { return deque.front().second; }
  [[nodiscard]] std::pair<int64_t, T> GetMinPair() const { return deque.front(); }

  void PopFront() {
    if (not deque.empty() and first == GetMinIndex()) {
      deque.pop_front();
    }
    ++first;
  }
};

class KmerMinimizerWindow {
 public:
  struct FreqHash {
    int64_t freq{0};
    Config::HashParams::htype hash{0};
    Config::HashParams::htype fhash{0};

    // bool operator==(const FreqHashPos &rhs) const { pos == rhs.pos; }
    // bool operator!=(const FreqHashPos &rhs) { return not operator==(rhs); }
    bool operator<(const FreqHash &rhs) const { return freq < rhs.freq or (freq == rhs.freq and hash < rhs.hash); }
    bool operator>(const FreqHash &rhs) const { return rhs.operator<(*this); }
  };

 private:
  MinQueue<FreqHash, std::greater<>> queue;

  static std::vector<FreqHash> GetInitWindow(KWH<Config::HashParams::htype> &kwh, const int64_t window_size,
                                             const kmer_index::KmerIndex::KmerCounter &counter) {
    std::vector<FreqHash> init_window;
    for (int i = 0; i < window_size; ++i, kwh = kwh.next()) {
      init_window.push_back({counter.at(kwh.hash()), kwh.hash(), kwh.get_fhash()});
      VERIFY(kwh.hasNext());
    }
    return init_window;
  }

  void PushBack(const int32_t freq, const Config::HashParams::htype hash, const Config::HashParams::htype fhash,
                const int64_t pos) {
    queue.PushBack(pos, {freq, hash, fhash});
  }
  void PopFront() { queue.PopFront(); }

 public:
  explicit KmerMinimizerWindow(const std::vector<FreqHash> &init_window) {
    for (auto it = init_window.cbegin(); it != init_window.cend(); ++it) {
      PushBack(it->freq, it->hash, it->fhash, it - init_window.cbegin());
    }
  }

  KmerMinimizerWindow(KWH<Config::HashParams::htype> &kwh, const int64_t window_size,
                      const kmer_index::KmerIndex::KmerCounter &counter)
      : KmerMinimizerWindow(GetInitWindow(kwh, window_size, counter)) {}

  void Add(const int32_t freq, const Config::HashParams::htype hash, const Config::HashParams::htype fhash,
           const int64_t pos) {
    PushBack(freq, hash, fhash, pos);
    PopFront();
  }

  [[nodiscard]] FreqHash GetMinimizer() const { return queue.GetMin(); }
  [[nodiscard]] int64_t GetMinimizerPos() const { return queue.GetMinIndex(); }
};
}// End namespace veritymap::kmer_index::kmer_window