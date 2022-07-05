//
// Created by Andrey Bzikadze on 06/14/22.
//

#pragma once

namespace veritymap::kmer_index {

class KmerIndex {
 public:
  using Kmer2PosSingle = std::unordered_map<Config::HashParams::htype, std::vector<int64_t>>;
  using Kmer2Pos = std::vector<Kmer2PosSingle>;
  using KmerCounter = std::unordered_map<Config::HashParams::htype, int64_t>;

 private:
  Kmer2Pos kmer2pos;
  KmerCounter counter;
  const std::vector<Contig> &ctgs;

 public:
  KmerIndex(Kmer2Pos kmer2pos, KmerCounter counter, const std::vector<Contig> &ctgs)
      : kmer2pos{std::move(kmer2pos)},
        counter{std::move(counter)},
        ctgs{ctgs} {
    VERIFY(this->ctgs.size() == this->kmer2pos.size());
  }

  KmerIndex(const std::vector<Contig> &ctgs, const std::filesystem::path &input_fn) : ctgs{ctgs} {
    std::unordered_map<std::string, int64_t> name2index;
    for (auto it = ctgs.begin(); it != ctgs.end(); ++it) { name2index.emplace(it->id, it - ctgs.begin()); }
    std::string name;
    int64_t pos;
    Config::HashParams::htype hash;
    kmer2pos.resize(name2index.size());
    std::ifstream is(input_fn);
    while (is >> name >> pos >> hash) {
      ++counter[hash];
      kmer2pos[name2index[name]][hash].push_back(pos);
    }
  }

  std::vector<int64_t> NSolidKmers() const {
    std::vector<int64_t> cnt;
    for (const auto &k2p : kmer2pos) { cnt.push_back(k2p.size()); }
    return cnt;
  }

  friend std::ostream &operator<<(std::ostream &os, const KmerIndex &index);

  Kmer2PosSingle &operator[](const int64_t i) { return kmer2pos[i]; }
  [[nodiscard]] const Kmer2PosSingle &operator[](const int64_t i) const { return kmer2pos[i]; }

  [[nodiscard]] int64_t GetCount(const Config::HashParams::htype &hash) const {
    auto it = counter.find(hash);
    return it == counter.end() ? 0 : it->second;
  }

  [[nodiscard]] int64_t GetCount(const Config::HashParams::htype &hash, const int64_t i) const {
    const Kmer2PosSingle &kmer2pos_single = kmer2pos.at(i);
    auto it = kmer2pos_single.find(hash);
    return it == kmer2pos_single.end() ? 0 : it->second.size();
  }

  [[nodiscard]] const std::vector<int64_t> *GetPos(const Config::HashParams::htype &hash, const int64_t i) const {
    auto it = kmer2pos.at(i).find(hash);
    return it != kmer2pos.at(i).end() ? &(it->second) : nullptr;
  }

  friend class HighFreqUniqueKmersFilterer;
};

std::ostream &operator<<(std::ostream &os, const KmerIndex &index) {
  for (auto it = index.kmer2pos.cbegin(); it != index.kmer2pos.cend(); ++it) {
    const Contig &contig = index.ctgs.at(it - index.kmer2pos.cbegin());
    for (const auto &[hash, pos] : *it) {
      for (const int64_t p : pos) { os << contig.id << "\t" << p << "\t" << hash << "\n"; }
    }
  }
  return os;
}

class HighFreqUniqueKmersFilterer {
  int64_t nthreads{1};
  const RollingHash<Config::HashParams::htype> &hasher;
  Config::CommonParams common_params;
  Config::KmerIndexerParams kmer_indexer_params;
  logging::Logger &logger;

 public:
  HighFreqUniqueKmersFilterer(const int64_t nthreads, const RollingHash<Config::HashParams::htype> &hasher,
                              const Config::CommonParams &common_params,
                              const Config::KmerIndexerParams &kmer_indexer_params, logging::Logger &logger)
      : nthreads{nthreads},
        hasher{hasher},
        common_params{common_params},
        kmer_indexer_params{kmer_indexer_params},
        logger{logger} {}

  void Filter(KmerIndex &kmer_index, const std::vector<Contig> &contigs) {
    // ban unique k-mers in assembly that have unusually high coverage

    logger.info() << "Counting unique k-mers from the target...\n";
    std::unordered_map<Config::HashParams::htype, std::atomic<size_t>> unique_kmers;
    for (const auto &[hash, cnt] : kmer_index.counter) {
      if (cnt == 1) {
        unique_kmers.emplace(hash, 0);
      }
    }
    logger.info() << "There are " << unique_kmers.size() << " unique k-mers in the target\n";

    std::function<void(const Contig &)> process_read = [&unique_kmers, this](const Contig &contig) {
      if (contig.size() < hasher.k) {
        return;
      }
      KWH<Config::HashParams::htype> kwh(hasher, contig.seq, 0);
      while (true) {
        if (!kwh.hasNext()) {
          break;
        }
        kwh = kwh.next();
        const Config::HashParams::htype fhash = kwh.get_fhash();
        const Config::HashParams::htype rhash = kwh.get_rhash();
        for (const Config::HashParams::htype hash : std::vector<Config::HashParams::htype>{fhash, rhash}) {
          auto kmer_it = unique_kmers.find(hash);
          if (kmer_it != unique_kmers.end()) {
            ++(kmer_it->second);
            break;
          }
        }
      }
    };
    process_in_parallel(contigs, process_read, nthreads, true);

    logger.info() << "Finished counting frequences of unique k-mers in the queries...\n";

    auto [mean, stddev] = [&unique_kmers]() -> std::pair<double, double> {
      auto it = unique_kmers.begin();
      double mean = it->second;
      double variance = 0;
      for (int k = 1; it != unique_kmers.end(); ++it, ++k) {
        double mean_pre = mean;
        const auto &val = it->second;
        mean += (val - mean) / k;
        variance += (val - mean) * (val - mean_pre);
      }
      return {mean, std::sqrt(variance / unique_kmers.size())};
    }();

    logger.info() << "Mean (std) multiplicity of a unique k-mer = " << mean << " (" << stddev << ")\n";

    const uint max_read_freq = mean + kmer_indexer_params.careful_upper_bnd_cov_mult * stddev;
    logger.info() << "Max solid k-mer frequency in reads " << max_read_freq << "\n";

    uint64_t n{0};
    for (auto &[hash, cnt] : unique_kmers) {
      if (cnt > max_read_freq) {
        kmer_index.counter.erase(hash);
        for (KmerIndex::Kmer2PosSingle &kmer2pos_single : kmer_index.kmer2pos) {
          auto it = kmer2pos_single.find(hash);
          if (it != kmer2pos_single.end()) {
            kmer2pos_single.erase(it);
            kmer_index.counter.erase(hash);
            break;
          }
        }
        ++n;
      }
    }
    logger.info() << "Filtered " << n << " high multiplicity k-mers\n";
  }
};

}// End namespace veritymap::kmer_index