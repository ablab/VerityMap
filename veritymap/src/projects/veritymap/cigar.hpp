//
// Created by Andrey Bzikadze on 03/02/21.
//

#pragma once

#include <sstream>
#include <string>
#include <vector>

#include "ksw2/ksw2.h"
#include "sequences/sequence.hpp"

namespace veritymap::cigar_utils {

enum class CigarMode { M,
                       I,
                       D,
                       S };

inline char cigar_mode2str(const CigarMode &fragment);

struct CigarFragment {
  size_t length{0};
  CigarMode mode{};
};

class Cigar {
  std::vector<CigarFragment> cigar_vec;

 public:
  explicit Cigar(const ksw_extz_t &ez);
  Cigar(const size_t length, const CigarMode mode) : cigar_vec{{length, mode}} {}

  Cigar() = default;
  Cigar(Cigar &) = default;
  Cigar(Cigar &&) = default;
  Cigar &operator=(const Cigar &) = default;
  Cigar &operator=(Cigar &&) = default;
  ~Cigar() = default;

  [[nodiscard]] bool empty() const { return cigar_vec.empty(); }

  [[nodiscard]] const std::vector<CigarFragment> &get_cigar_vec() const;

  void extend(size_t length, CigarMode mode);

  void extend(Cigar cigar);

  [[nodiscard]] size_t query_length() const;

  [[nodiscard]] size_t target_length() const;

  [[nodiscard]] int nmismatches(const Sequence &target, const Sequence &query) const;

  [[nodiscard]] int alignment_length() const;

  [[nodiscard]] double identity(const Sequence &target, const Sequence &query) const;

  std::pair<size_t, size_t> trim(const CigarMode &mode) {
    size_t left_trim = {0}, right_trim{0};

    if (cigar_vec.empty()) { return {left_trim, right_trim}; }
    if (cigar_vec.front().mode == mode) {
      left_trim = cigar_vec.front().length;
      cigar_vec.erase(cigar_vec.begin());
    }

    if (cigar_vec.empty()) { return {left_trim, right_trim}; }
    if (cigar_vec.back().mode == mode) {
      right_trim = cigar_vec.back().length;
      cigar_vec.pop_back();
    }
    return {left_trim, right_trim};
  }

  void soft_clip() {
    if (cigar_vec.empty()) { return; }
    if (cigar_vec.front().mode == CigarMode::I) { cigar_vec.front().mode = CigarMode::S; }
    if (cigar_vec.back().mode == CigarMode::I) { cigar_vec.back().mode = CigarMode::S; }
  }
};

std::ostream &operator<<(std::ostream &os, const Cigar &cigar);

}// namespace veritymap::cigar_utils