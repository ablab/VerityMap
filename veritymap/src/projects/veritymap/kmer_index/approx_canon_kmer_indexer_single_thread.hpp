//
// Created by Andrey Bzikadze on 06/02/22.
//

#pragma once

namespace veritymap::kmer_index::approx_canon_single_thread_kmer_indexer {

template<typename htype>
class ApproxCanonSingleThreadKmerIndexer {
  const RollingHash<htype> &hasher;
  Config::CommonParams common_params;
  Config::KmerIndexerParams kmer_indexer_params;

 public:
  ApproxCanonSingleThreadKmerIndexer(const RollingHash<htype> &hasher, const Config::CommonParams &common_params,
                                     const Config::KmerIndexerParams &kmer_indexer_params)
      : hasher{hasher},
        common_params{common_params},
        kmer_indexer_params{kmer_indexer_params} {}

  [[nodiscard]] KmerIndexes extract(const std::vector<Contig> &contigs, const std::vector<Contig> &readset,
                                    logging::Logger &logger) const {
    int64_t tot_len{0};
    for (const Contig &contig : contigs) { tot_len += contig.size(); }
    const cms_utils::CMSParams kCmsParams(common_params, kmer_indexer_params, tot_len, 1);
    sketch::cm::ccm_t cms(kCmsParams.nbits, kCmsParams.l2sz, kCmsParams.nhash);

    for (const Contig &contig : contigs) {
      if (contig.size() < hasher.k) {
        continue;
      }
      KWH<Config::HashParams::htype> kwh(hasher, contig.seq, 0);
      while (true) {
        htype hash = kwh.hash();
        if (cms.est_count(hash) <= kmer_indexer_params.max_rare_cnt_target) {
          cms.add(kwh.hash());
        }

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
        if (cms.est_count(kwh.hash()) <= kmer_indexer_params.max_rare_cnt_target) {
          kmer_index[kwh.get_fhash()].emplace_back(kwh.pos);
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

}// namespace veritymap::kmer_index::approx_canon_single_thread_kmer_indexer