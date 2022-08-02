//
// Created by Andrey Bzikadze on 06/15/22.
//

#pragma once

#include "kmer_index.hpp"

namespace veritymap::indexed_contigs {

class IndexedContigs {
  const std::vector<Contig>& contigs;
  const RollingHash<Config::HashParams::htype>& hasher;
  kmer_index::KmerIndex index;

  std::vector<std::pair<int64_t, int64_t>> GetNoSolidRegionsPerContig(const int i, const int max_dist,
                                                                      const int64_t k) const {
    const Contig& contig = contigs.at(i);
    std::vector<int64_t> pos;
    for (const auto& [hash, kmer_pos] : index[i]) {
      for (const int64_t p : kmer_pos) { pos.emplace_back(p); }
    }
    std::sort(pos.begin(), pos.end());

    int64_t prev_p{k};
    std::vector<std::pair<int64_t, int64_t>> norare_regions;
    for (const auto& p : pos) {
      if (p <= prev_p) {
        continue;
      }
      if (p - prev_p > max_dist) {
        norare_regions.emplace_back(prev_p, p);
      }
      prev_p = p + k;
    }
    return norare_regions;
  }

 public:
  IndexedContigs(const std::vector<Contig>& contigs, const RollingHash<Config::HashParams::htype>& hasher,
                 kmer_index::KmerIndex index)
      : contigs{contigs},
        hasher{hasher},
        index{std::move(index)} {}

  void Summary(logging::Logger& logger) const {
    std::vector<int64_t> n_solid_kmers = index.NSolidKmers();
    for (auto it = contigs.cbegin(); it != contigs.cend(); ++it) {
      logger.info() << "Sequence " << it->id << ", # Solid kmers = " << n_solid_kmers[it - contigs.cbegin()]
                    << std::endl;
    }
  }

  std::ostream& NoSolidRegions2Bed(const int64_t max_dist, const int64_t k, std::ostream& os) const {
    for (int64_t i = 0; i < contigs.size(); ++i) {
      const std::vector<std::pair<int64_t, int64_t>> norare_regions = GetNoSolidRegionsPerContig(i, max_dist, k);
      const Contig& contig = contigs[i];
      for (const auto [s, e] : norare_regions) {
        os << contig.id << "\t" << s << "\t" << e << "\t" << e - s << " bp\n";
      }
    }
    return os;
  }

  [[nodiscard]] const std::vector<Contig>& Contigs() const { return contigs; }
  [[nodiscard]] int64_t Size() const { return contigs.size(); }

  [[nodiscard]] const kmer_index::KmerIndex& Index() const { return index; }
};

}// namespace veritymap::indexed_contigs