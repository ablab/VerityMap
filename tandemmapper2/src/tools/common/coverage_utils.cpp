//
// Created by Andrey Bzikadze on 10/19/21.
//

#include "coverage_utils.hpp"

using namespace tools::common::coverage_utils;

double tools::common::coverage_utils::get_coverage(const std::vector<Contig>& contigs,
                                                   const std::vector<Contig>& readset) {
  uint64_t cnt_len{0};
  for (const Contig& contig : contigs) {
    cnt_len += contig.size();
  }

  uint64_t reads_len{0};
  for (const Contig& read : readset) {
    reads_len += read.size();
  }
  return static_cast<double>(reads_len) / static_cast<double>(cnt_len);
}
