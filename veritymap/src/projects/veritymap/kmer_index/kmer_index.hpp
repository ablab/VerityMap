//
// Created by Andrey Bzikadze on 2/22/21.
//

#pragma once

#include <unordered_set>

#include "../config/config.hpp"
#include "../hash_utils.hpp"
#include "../rolling_hash.hpp"

namespace veritymap::kmer_index {

using KmerIndex = std::unordered_map<Config::HashParams::htype, std::vector<size_t>>;
using KmerIndexes = std::vector<KmerIndex>;

using Counter = std::unordered_map<Config::HashParams::htype, size_t>;

using Counters = std::vector<Counter>;

Counters _get_counters(const RollingHash<Config::HashParams::htype>& hasher,
                       const std::vector<Sequence>& seqs) {
  Counters counters;
  for (auto it = seqs.cbegin(); it != seqs.cend(); ++it) {
    const Sequence& seq = *it;
    std::unordered_map<Config::HashParams::htype, size_t>& counter{counters.emplace_back()};
    if (seq.size() < hasher.k) {
      continue;
    }
    const size_t i = it - seqs.cbegin();
    KWH<Config::HashParams::htype> kwh(hasher, seq, 0);
    while (true) {
      counter[kwh.get_fhash()] += 1;
      if (!kwh.hasNext()) {
        break;
      }
      kwh = kwh.next();
    }
  }
  return counters;
}

Counter _get_readset_counter(const RollingHash<Config::HashParams::htype>& hasher,
                             const std::vector<Contig>& contigs) {
  std::vector<Sequence> seqs;
  for (const auto& contig : contigs) {
    seqs.emplace_back(contig.seq);
  }
  Counter counter;
  for (auto it = seqs.cbegin(); it != seqs.cend(); ++it) {
    const Sequence& seq = *it;
    if (seq.size() < hasher.k) {
      continue;
    }
    KWH<Config::HashParams::htype> kwh(hasher, seq, 0);
    while (true) {
      counter[kwh.hash()] += 1;
      if (!kwh.hasNext()) {
        break;
      }
      kwh = kwh.next();
    }
  }
  return counter;
}

std::unordered_map<Config::HashParams::htype, std::unordered_set<size_t>>
_get_hash2seqs(const Counters& counters) {
  std::unordered_map<Config::HashParams::htype, std::unordered_set<size_t>> hash2seqs;
  for (auto it = counters.cbegin(); it != counters.cend(); ++it) {
    for (const auto& [hash, cnt] : *it) {
      hash2seqs[hash].insert(it - counters.cbegin());
    }
  }
  return hash2seqs;
}

KmerIndexes _get_rare_kmers_from_counter(const RollingHash<Config::HashParams::htype>& hasher,
                                         const std::vector<Sequence>& seqs,
                                         const std::unordered_map<
                                             Config::HashParams::htype,
                                             std::unordered_set<size_t>>& hash2seqs,
                                         const Counters& counters,
                                         const size_t max_rare_cnt) {
  KmerIndexes rare_kmers_indexes;
  for (auto it = seqs.cbegin(); it != seqs.cend(); ++it) {
    KmerIndex& rare_kmers_index{rare_kmers_indexes.emplace_back()};
    const Sequence& seq = *it;
    const Counter& counter = counters.at(it - seqs.cbegin());
    if (seq.size() < hasher.k) {
      continue;
    }
    KWH<Config::HashParams::htype> kwh(hasher, seq, 0);
    while (true) {
      Config::HashParams::htype fhash{kwh.get_fhash()};
      Config::HashParams::htype rhash{kwh.get_rhash()};
      if ((hash2seqs.at(fhash).size() == 1) and (not hash2seqs.contains(rhash)) and (counter.at(fhash) <= max_rare_cnt)) {
        rare_kmers_index[fhash].emplace_back(kwh.pos);
      }
      if (!kwh.hasNext()) {
        break;
      }
      kwh = kwh.next();
    }
  }
  return rare_kmers_indexes;
}

KmerIndexes get_rare_kmers(const std::vector<Sequence>& sequences,
                           const RollingHash<Config::HashParams::htype>& hasher,
                           const size_t max_rare_cnt) {
  const Counters counters = _get_counters(hasher, sequences);
  const std::unordered_map<Config::HashParams::htype, std::unordered_set<size_t>> hash2seqs =
      _get_hash2seqs(counters);
  KmerIndexes rare_kmers_indexes = _get_rare_kmers_from_counter(hasher,
                                                                sequences,
                                                                hash2seqs,
                                                                counters,
                                                                max_rare_cnt);
  return rare_kmers_indexes;
}

KmerIndexes get_rare_kmers(const std::vector<Contig>& contigs,
                           const RollingHash<Config::HashParams::htype>& hasher,
                           const size_t max_rare_cnt) {
  std::vector<Sequence> seqs;
  for (const auto& contig : contigs) {
    seqs.emplace_back(contig.seq);
  }
  return get_rare_kmers(seqs, hasher, max_rare_cnt);
}

KmerIndex get_rare_kmers(const Sequence& sequence,
                         const RollingHash<Config::HashParams::htype>& hasher,
                         const size_t max_rare_cnt) {
  const std::vector<Sequence> sequences{sequence};
  return get_rare_kmers(sequences, hasher, max_rare_cnt).front();
}

std::vector<std::pair<size_t, size_t>> _get_norare_regions(const KmerIndex& kmer_index,
                                                           const size_t max_dist, const size_t k,
                                                           const size_t seq_len) {
  VERIFY(seq_len > 0);
  std::vector<size_t> pos{seq_len - 1};
  for (const auto& [hash, kmer_pos] : kmer_index) {
    for (const size_t p : kmer_pos) {
      pos.emplace_back(p);
    }
  }
  std::sort(pos.begin(), pos.end());

  size_t prev_p{k};
  std::vector<std::pair<size_t, size_t>> norare_regions;
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

class IndexedContig {
  const Contig& contig;
  const RollingHash<Config::HashParams::htype>& hasher;
  size_t max_rare_cnt;
  KmerIndex kmer_index;

 public:
  IndexedContig(const Contig& contig_,
                const RollingHash<Config::HashParams::htype>& hasher_,
                const size_t max_rare_cnt_,
                KmerIndex kmer_index = {}) : contig{contig_},
                                             hasher{hasher_},
                                             max_rare_cnt{max_rare_cnt_},
                                             kmer_index{std::move(kmer_index)} {}

  [[nodiscard]] const Contig& get_contig() const { return contig; }
  [[nodiscard]] size_t get_contig_size() const { return contig.size(); }
  const RollingHash<Config::HashParams::htype>& get_hasher() const { return hasher; }
  [[nodiscard]] size_t get_max_rare_cnt() const { return max_rare_cnt; }
  const KmerIndex& get_kmer_index() const { return kmer_index; }

  KmerIndex& get_kmer_index() { return kmer_index; }

  [[nodiscard]] std::vector<std::pair<size_t, size_t>>
  get_norare_regions(const size_t max_dist, const size_t k) const {
    return _get_norare_regions(kmer_index, max_dist, k, get_contig_size());
  }

  std::ostream& norare_regions2bam(const size_t max_dist, const size_t k, std::ostream& os) const {
    const std::vector<std::pair<size_t, size_t>> norare_regions = get_norare_regions(max_dist, k);
    for (const auto [s, e] : norare_regions) {
      os << contig.id << "\t" << s << "\t" << e << "\t" << e - s << " bp\n";
    }
    return os;
  }

  std::ostream& export_index(std::ostream& os) const {
    for (const auto& [hash, pos] : kmer_index) {
      for (const size_t p : pos) {
        os << contig.id << "\t" << p << "\t" << hash << "\n";
      }
    }
    return os;
  }
};

using IndexedContigs = std::vector<IndexedContig>;

std::ostream& export_index(std::ostream& os, const IndexedContigs& kmer_indexes) {
  for (const IndexedContig& index : kmer_indexes) {
    index.export_index(os);
  }
  return os;
}

IndexedContigs import_index(const std::vector<Contig>& contigs,
                            const RollingHash<Config::HashParams::htype>& hasher,
                            const size_t max_rare_cnt,
                            std::istream& is) {
  IndexedContigs icontigs;
  std::unordered_map<std::string, size_t> name2ctg;
  for (auto it = contigs.begin(); it != contigs.end(); ++it) {
    icontigs.emplace_back(*it, hasher, max_rare_cnt);
    name2ctg[it->id] = it - contigs.begin();
  }

  std::string name;
  size_t pos;
  Config::HashParams::htype hash;
  while (is >> name >> pos >> hash) {
    icontigs[name2ctg[name]].get_kmer_index()[hash].emplace_back(pos);
  }
  return icontigs;
}

IndexedContigs get_exact_indexed_contigs(const std::vector<Contig>& contigs,
                                         const RollingHash<Config::HashParams::htype>& hasher,
                                         size_t max_rare_cnt) {
  IndexedContigs indexed_contigs;
  KmerIndexes rare_kmers_indexes = get_rare_kmers(contigs, hasher, max_rare_cnt);
  for (auto it = rare_kmers_indexes.begin(); it != rare_kmers_indexes.end(); ++it) {
    const Contig& contig = contigs.at(it - rare_kmers_indexes.begin());
    indexed_contigs.emplace_back(contig, hasher, max_rare_cnt, std::move(*it));
  }
  return indexed_contigs;
}

std::ostream& norare_regions2bam(const IndexedContigs& indexed_contigs,
                                 const size_t max_dist, const size_t k, std::ostream& os) {
  for (const IndexedContig& indexed_contig : indexed_contigs) {
    indexed_contig.norare_regions2bam(max_dist, k, os);
  }
  return os;
}

}// End namespace veritymap::kmer_index