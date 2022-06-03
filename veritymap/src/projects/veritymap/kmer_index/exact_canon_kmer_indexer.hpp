//
// Created by Andrey Bzikadze on 06/02/22.
//

#pragma once

namespace veritymap::kmer_index::exact_canon_kmer_indexer {

template<typename htype>
class ExactCanonKmerIndexer {
  const RollingHash<htype> &hasher;
  Config::CommonParams common_params;
  Config::KmerIndexerParams kmer_indexer_params;

 public:
  ExactCanonKmerIndexer(const RollingHash<htype> &hasher, const Config::CommonParams &common_params,
                        const Config::KmerIndexerParams &kmer_indexer_params)
      : hasher{hasher},
        common_params{common_params},
        kmer_indexer_params{kmer_indexer_params} {}

  [[nodiscard]] KmerIndexes extract(const std::vector<Contig> &contigs, const std::vector<Contig> &readset,
                                    logging::Logger &logger) const {
    Counter counter;

    for (const Contig &contig : contigs) {
      if (contig.size() < hasher.k) {
        continue;
      }
      KWH<Config::HashParams::htype> kwh(hasher, contig.seq, 0);
      while (true) {
        counter[kwh.hash()] += 1;
        if (!kwh.hasNext()) {
          break;
        }
        kwh = kwh.next();
      }
    }

    KmerIndexes kmer_indexes;
    for (const Contig &contig : contigs) {
      KmerIndex &kmer_index = kmer_indexes.emplace_back();
      if (contig.size() < hasher.k) {
        continue;
      }
      KWH<Config::HashParams::htype> kwh(hasher, contig.seq, 0);
      while (true) {
        if (counter.at(kwh.hash()) <= kmer_indexer_params.max_rare_cnt_target) {
          kmer_index[kwh.get_fhash()].template emplace_back(kwh.pos);
        }
        if (!kwh.hasNext()) {
          break;
        }
        kwh = kwh.next();
      }
    }

    // TODO allow banning high frequency unique kmers
    // logger.info() << "Filtering high multiplicity unique k-mers\n";
    // BanHighFreqUniqueKmers(contigs, readset, kmer_indexes, logger);

    return kmer_indexes;
  }
};

}// namespace veritymap::kmer_index::exact_canon_kmer_indexer