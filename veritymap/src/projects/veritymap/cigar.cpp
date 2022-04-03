//
// Created by Andrey Bzikadze on 03/02/21.
//

#include "cigar.hpp"

#include "sequences/verify.hpp"

using namespace veritymap::cigar_utils;

inline char veritymap::cigar_utils::cigar_mode2str(const CigarMode& fragment) {
  if (fragment == CigarMode::M) {
    return 'M';
  } else if (fragment == CigarMode::I) {
    return 'I';
  } else if (fragment == CigarMode::D) {
    return 'D';
  } else {
    VERIFY(fragment == CigarMode::S);
    return 'S';
  }
}

Cigar::Cigar(const ksw_extz_t& ez) {
  for (size_t i = 0; i < ez.n_cigar; ++i) {
    const auto mode = static_cast<CigarMode>(ez.cigar[i] & 0xf);
    const size_t length{ez.cigar[i] >> 4};
    cigar_vec.push_back({length, mode});
  }
}

const std::vector<CigarFragment>& Cigar::get_cigar_vec() const { return cigar_vec; }

void Cigar::extend(const size_t length, const CigarMode mode) {
  if (empty() or cigar_vec.back().mode != mode) {
    cigar_vec.push_back({length, mode});
  } else {
    cigar_vec.back().length += length;
  }
}

void Cigar::extend(Cigar cigar) {
  if ((not empty()) and (not cigar.empty()) and (cigar_vec.back().mode == cigar.cigar_vec.front().mode)) {
    cigar.cigar_vec.front().length += cigar_vec.back().length;
    cigar_vec.pop_back();
  }
  cigar_vec.insert(cigar_vec.end(), std::make_move_iterator(cigar.cigar_vec.begin()),
                   std::make_move_iterator(cigar.cigar_vec.end()));
}

[[nodiscard]] size_t Cigar::query_length() const {
  size_t length{0};
  for (const CigarFragment& fragment : cigar_vec) {
    if (fragment.mode != CigarMode::D) {
      length += fragment.length;
    }
  }
  return length;
}

[[nodiscard]] size_t Cigar::target_length() const {
  size_t length{0};
  for (const CigarFragment& fragment : cigar_vec) {
    if (fragment.mode != CigarMode::I and fragment.mode != CigarMode::S) {
      length += fragment.length;
    }
  }
  return length;
}

[[nodiscard]] int Cigar::nmismatches(const Sequence& target, const Sequence& query) const {
  int nmism{0};
  size_t t{0}, q{0};
  for (const CigarFragment& fragment : cigar_vec) {
    if (fragment.mode == CigarMode::I or fragment.mode == CigarMode::S) {
      q += fragment.length;
      nmism += fragment.length;
    } else if (fragment.mode == CigarMode::D) {
      t += fragment.length;
      nmism += fragment.length;
    } else {
      VERIFY(fragment.mode == CigarMode::M);
      for (size_t i = 0; i < fragment.length; ++i) {
        if (target[t] != query[q]) {
          ++nmism;
        }
        ++t, ++q;
      }
    }
  }
  return nmism;
}

int Cigar::alignment_length() const {
  int length{0};
  for (const CigarFragment& fragment : cigar_vec) { length += fragment.length; }
  return length;
}

double Cigar::identity(const Sequence& target, const Sequence& query) const {
  const int al_len = alignment_length();
  if (al_len > 0) {
    const double identity = 1. - static_cast<double>(nmismatches(target, query)) / al_len;
    VERIFY(identity >= 0);
    return identity;
  }
  return 1.;
}

std::ostream& veritymap::cigar_utils::operator<<(std::ostream& os, const Cigar& cigar) {
  for (const CigarFragment& fragment : cigar.get_cigar_vec()) {
    os << fragment.length << cigar_mode2str(fragment.mode);
  }
  return os;
}