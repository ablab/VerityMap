//
// Created by Andrey Bzikadze on 03/15/21.
//

#pragma once

#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>

#include "../config/config.hpp"
#include "../rolling_hash.hpp"
#include "approx_canon_kmer_indexer.hpp"
#include "approx_kmer_indexer.hpp"
#include "bloom/bloom.hpp"
#include "include/sketch/ccm.h"

namespace veritymap::kmer_index {

using Counter = std::unordered_map<Config::HashParams::htype, size_t>;

kmer_index::IndexedContigs get_indexed_targets(const std::vector<Contig> &queries, const std::vector<Contig> &targets,
                                               const std::filesystem::path &outdir,
                                               const RollingHash<Config::HashParams::htype> &hasher,
                                               const size_t nthreads, logging::Logger &logger,
                                               const std::filesystem::path &index_path,
                                               const Config::CommonParams &common_params,
                                               const Config::KmerIndexerParams &kmer_indexer_params) {
  using htype = Config::HashParams::htype;
  const auto kmer_indexes_fn = outdir / "kmer_indexes.tsv";

  if (not index_path.empty()) {
    logger.info() << "Importing kmer index from " << index_path << std::endl;
    std::ifstream kmer_indexes_is(index_path);
    IndexedContigs indexed_targets =
        import_index(targets, hasher, kmer_indexer_params.max_rare_cnt_target, kmer_indexes_is);
    kmer_indexes_is.close();
    logger.info() << "Finished importing kmer index from " << index_path << std::endl;
    std::filesystem::copy_file(index_path, kmer_indexes_fn, std::filesystem::copy_options::overwrite_existing);
    return indexed_targets;
  }

  std::vector<KmerIndex> kmers_indexes = [&queries, &targets, &nthreads, &hasher, &common_params, &kmer_indexer_params,
                                          &index_path, &logger] {
    if (kmer_indexer_params.strategy == Config::KmerIndexerParams::Strategy::exact) {
      logger.info() << "Getting exact kmer indexes..." << std::endl;
      std::vector<KmerIndex> kmers_indexes = get_rare_kmers(targets, hasher, kmer_indexer_params.max_rare_cnt_target);
      logger.info() << "Finished getting exact kmer indexes" << std::endl;
      return kmers_indexes;
    }
    if (kmer_indexer_params.strategy == Config::KmerIndexerParams::Strategy::approximate) {
      logger.info() << "Getting approximate kmer indexes..." << std::endl;

      const approx_kmer_indexer::ApproxKmerIndexer kmer_indexer(nthreads, hasher, common_params, kmer_indexer_params);
      KmerIndexes kmers_indexes = kmer_indexer.extract(targets, queries, logger);
      logger.info() << "Finished getting approximate kmer indexes" << std::endl;
      return kmers_indexes;
    } else {
      VERIFY(kmer_indexer_params.strategy == Config::KmerIndexerParams::Strategy::approximate_canon);
      logger.info() << "Getting approximate kmer indexes (canonical variant)..." << std::endl;
      const approx_canon_kmer_indexer::ApproxCanonKmerIndexer kmer_indexer(nthreads, hasher, common_params,
                                                                           kmer_indexer_params);
      KmerIndexes kmers_indexes = kmer_indexer.extract(targets, queries, logger);
      logger.info() << "Finished getting approximate (canonical variant) kmer indexes" << std::endl;
      return kmers_indexes;
    }
  }();

  IndexedContigs indexed_targets;
  for (auto it = kmers_indexes.begin(); it != kmers_indexes.end(); ++it) {
    const Contig &target = targets.at(it - kmers_indexes.begin());
    indexed_targets.emplace_back(target, hasher, kmer_indexer_params.max_rare_cnt_target, std::move(*it));
  }

  std::ofstream kmer_indexes_os(kmer_indexes_fn);
  veritymap::kmer_index::export_index(kmer_indexes_os, indexed_targets);
  kmer_indexes_os.close();
  logger.info() << "Kmer indexes are exported to " << kmer_indexes_fn << std::endl;

  for (const auto &indexed_target : indexed_targets) {
    logger.info() << "Target " << indexed_target.get_contig().id
                  << ", # Rare kmers = " << indexed_target.get_kmer_index().size() << std::endl;
  }

  return indexed_targets;
}
}// End namespace veritymap::kmer_index