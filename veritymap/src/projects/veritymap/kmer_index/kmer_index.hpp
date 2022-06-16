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

}// End namespace veritymap::kmer_index