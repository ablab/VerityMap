//
// Created by Andrey Bzikadze on 2/19/21.
//

#pragma once

#include <common/parallel.h>

#include <common/logging.hpp>
#include <mutex>
#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>

#include "chaining.hpp"
#include "config/config.hpp"
#include "dp_scoring.hpp"
#include "hash_utils.hpp"
#include "kmer_index/kmer_index.hpp"
#include "kmer_index/target_indexer.hpp"
#include "matches.hpp"

namespace veritymap {

using Counter = std::unordered_map<Config::HashParams::htype, size_t>;

chaining::Chains _map_single_strand(const kmer_index::IndexedContig &indexed_target,
                                    const Contig &query,
                                    const dna_strand::Strand &query_strand,
                                    const Config &config) {
  const matches::Matches matches =
      matches::get_matches(indexed_target.get_contig(),
                           indexed_target.get_kmer_index(),
                           query,
                           query_strand,
                           indexed_target.get_hasher(),
                           config.kmer_indexer_params);
  if (matches.size() < config.chaining_params.min_matches) {
    return {};
  }

  const auto [scores, backtracks] = scoring::get_scores(matches, config.common_params, config.chaining_params);

  chaining::Chains chains = chaining::get_chains(indexed_target.get_contig(),
                                                 query, query_strand,
                                                 matches, scores, backtracks,
                                                 config.common_params, config.chaining_params);
  return chains;
}

std::optional<chaining::Chain>
_map_single(const Contig &query, const kmer_index::IndexedContigs &indexed_targets, const Config &config) {
  using score_type = typename Config::ChainingParams::score_type;

  chaining::Chains chains;
  for (const kmer_index::IndexedContig &indexed_target : indexed_targets) {
    chaining::Chains chains_f =
        _map_single_strand(indexed_target, query, dna_strand::Strand::forward, config);
    chaining::Chains chains_r =
        _map_single_strand(indexed_target, query, dna_strand::Strand::reverse, config);

    for (size_t i = 0; i < 2; ++i) {
      if (chains_f.size() > i) {
        chains.emplace_back(std::move(chains_f[i]));
      }
      if (chains_r.size() > i) {
        chains.emplace_back(std::move(chains_r[i]));
      }
    }
  }
  if (chains.empty())
    return std::nullopt;
  if (chains.size() == 1)
    return std::move(chains.front());

  auto pr_it = std::max_element(chains.begin(), chains.end(),
                                [](const auto &lhs, const auto &rhs) { return lhs.score < rhs.score; });
  Config::ChainingParams::score_type sc_score = 0;
  for (auto it = chains.begin(); it != chains.end(); ++it) {
    if (it == pr_it)
      continue;
    if (it->score > sc_score)
      sc_score = it->score;
  }
  const double top_score_prop = static_cast<double>(sc_score) / pr_it->score;
  if (top_score_prop < config.chaining_params.max_top_score_prop)
    return std::move(*pr_it);
  return std::nullopt;
}

void parallel_run(const kmer_index::IndexedContigs &indexed_targets,
                  std::vector<Contig> &queries,
                  const RollingHash<typename Config::HashParams::htype> &hasher,
                  const size_t nthreads,
                  const std::filesystem::path &chains_fn,
                  const std::filesystem::path &sam_fn,
                  const std::string &cmd,
                  const Config &config,
                  logging::Logger &logger) {
  std::mutex chainsMutex;
  using TargetQuery = std::tuple<const kmer_index::IndexedContig *, const Contig *>;

  std::ofstream chains_os(chains_fn);
  std::ofstream sam_os(sam_fn);
  for (const kmer_index::IndexedContig &itarget : indexed_targets) {
    const Contig &target = itarget.get_contig();
    sam_os << "@SQ\tSN:" << target.id << "\tLN:" << target.seq.size() << "\n";
  }
  sam_os << "@PG\tID:VerityMap\tPN:VerityMap\tVN:2.0\tCL:" << cmd << "\n";

  std::function<void(const Contig &query)> align_read =
      [&indexed_targets, &chainsMutex, &chains_os, &sam_os, &config](const Contig &query) {
        std::optional<chaining::Chain> chain = _map_single(query, indexed_targets, config);

        if (chain.has_value()) {
          std::string sam_record = chaining::chain2samrecord(chain.value(),
                                                             config.common_params,
                                                             config.chain2sam_params);
          chainsMutex.lock();
          chains_os << chain.value();
          sam_os << sam_record << "\n";
          chainsMutex.unlock();
        }
      };

  process_in_parallel(queries, align_read, nthreads, true);

  chains_os.close();
  sam_os.close();
}

void map(const std::filesystem::path &target_path,
         const std::filesystem::path &queries_path,
         const std::filesystem::path &outdir,
         const bool to_compress,
         const bool only_index,
         const bool careful_mode,
         const size_t nthreads,
         logging::Logger &logger,
         const std::string &cmd,
         const std::filesystem::path &index_path,
         Config config) {
    if (to_compress) {
        // TODO change that to a parameter call
        StringContig::needs_compressing = true;
    }

    std::optional<std::vector<Contig>> queries_optional;
    if (careful_mode) {
        queries_optional = io::SeqReader(queries_path).readAllContigs();
    }

    const RollingHash<Config::HashParams::htype>
        hasher(config.common_params.k, config.hash_params.base);

    io::SeqReader target_reader(target_path);
    std::vector<Contig> targets{target_reader.readAllContigs()};

    for (const Contig &target : targets) {
        logger.info() << "Target length " << target.seq.size() << ", name "
                      << target.id
                      << " from " << target_path << std::endl;
    }

  //const Counter query_counter = kmer_index::_get_readset_counter(hasher, queries);
  //logger.info() << query_counter.size();
  const kmer_index::IndexedContigs indexed_targets =
      kmer_index::get_indexed_targets(queries_optional,
                                      targets,
                                      outdir,
                                      hasher,
                                      nthreads,
                                      logger,
                                      index_path,
                                      config.common_params,
                                      config.kmer_indexer_params);

    const auto uncovered_fn = outdir/"norarekmers.bed";
    std::ofstream uncovered_os(uncovered_fn);
    kmer_index::norare_regions2bam(indexed_targets,
                                   config.kmer_indexer_params.min_uncovered_len,
                                   config.common_params.k,
                                   uncovered_os);
    logger.info() << "Finished exporting long (>= "
                  << config.kmer_indexer_params.min_uncovered_len << " bp) "
                  << "regions without rare k-mers to " << uncovered_fn
                  << std::endl;
    uncovered_os.close();

    if (only_index) {
        if (to_compress) {
            // TODO change that to a parameter call
            StringContig::needs_compressing = false;
        }
        return;
    }

    std::vector<Contig> queries = careful_mode ?
                                  std::move(queries_optional.value()) :
                                  io::SeqReader(queries_path).readAllContigs();

    logger.info() << "Queries from " << queries_path << ", total "
                  << queries.size() << " sequences " << std::endl;

    const auto chains_fn = outdir/"chains.tsv";
    const auto sam_fn = outdir/"alignments.sam";

    logger.info() << "Computing chains and sam records..." << std::endl;
    parallel_run(indexed_targets,
                 queries,
                 hasher,
                 nthreads,
                 chains_fn,
                 sam_fn,
                 cmd,
                 config,
                 logger);

    logger.info() << "Finished outputting chains to " << chains_fn
                  << " and sam records to " << sam_fn << std::endl;

  if (to_compress) {
    // TODO change that to a parameter call
    StringContig::needs_compressing = false;
  }
}
}// End namespace veritymap