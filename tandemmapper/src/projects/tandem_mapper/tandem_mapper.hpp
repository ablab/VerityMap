//
// Created by Andrey Bzikadze on 2/19/21.
//

#pragma once

#include <mutex>
#include <common/logging.hpp>
#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>
#include <common/parallel.h>

#include "hash_utils.hpp"
#include "kmer_index/kmer_index.hpp"
#include "kmer_index/target_indexer.hpp"
#include "matches.hpp"
#include "dp_scoring.hpp"
#include "chaining.hpp"
#include "config/config.hpp"


namespace tandem_mapper {

    using Counter = std::unordered_map<Config::HashParams::htype, size_t>;

    chaining::Chains _map_single_strand(const kmer_index::IndexedContig & indexed_target,
                                        const Contig & query,
                                        const dna_strand::Strand & query_strand,
                                        const Config & config) {
        const matches::Matches matches =
                matches::get_matches(indexed_target.get_contig(),
                                     indexed_target.get_kmer_index(),
                                     query,
                                     query_strand,
                                     indexed_target.get_hasher(),
                                     config.kmer_indexer_params.max_rare_cnt_query);
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
    _map_single(const kmer_index::IndexedContig & indexed_target, const Contig & query, const Config & config) {
        using score_type = typename Config::ChainingParams::score_type;

        chaining::Chains chains_f =
                _map_single_strand(indexed_target, query, dna_strand::Strand::forward, config);
        chaining::Chains chains_r =
                _map_single_strand(indexed_target, query, dna_strand::Strand::reverse, config);

        chaining::Chains chains;
        for (size_t i = 0; i < 2; ++i) {
            if (chains_f.size() > i) {
                chains.emplace_back(chains_f[i]);
            }
            if (chains_r.size() > i) {
                chains.emplace_back(chains_r[i]);
            }
        }
        if (chains.empty())
            return std::nullopt;
        if (chains.size() == 1)
            return chains.front();

        auto pr_it = std::max_element(chains.begin(), chains.end(),
                                      [](const auto & lhs, const auto & rhs) { return lhs.score < rhs.score; });
        Config::ChainingParams::score_type sc_score = 0;
        for (auto it = chains.begin(); it != chains.end(); ++it) {
            if (it == pr_it)
                continue;
            if (it->score > sc_score)
                sc_score = it->score;
        }
        const double top_score_prop = static_cast<double>(sc_score) / pr_it->score;
        if (top_score_prop < config.chaining_params.max_top_score_prop)
            return *pr_it;
        return std::nullopt;
    }

    void parallel_run(const kmer_index::IndexedContigs & indexed_targets,
                      std::vector<Contig> & queries,
                      const RollingHash<typename Config::HashParams::htype> & hasher,
                      const size_t nthreads,
                      const std::filesystem::path & chains_fn,
                      const std::filesystem::path & sam_fn,
                      const std::string & cmd,
                      const Config & config) {
        std::mutex chainsMutex;
        using TargetQuery = std::tuple<const kmer_index::IndexedContig *, const Contig *>;

        std::ofstream chains_os(chains_fn);
        std::ofstream sam_os(sam_fn);
        for (const kmer_index::IndexedContig & itarget : indexed_targets) {
            const Contig & target = itarget.get_contig();
            sam_os << "@SQ\tSN:" << target.id << "\tLN:" << target.seq.size() << "\n";
        }
        sam_os << "@PG\tID:tandemMapper2\tPN:tandemMapper2\tVN:2.0\tCL:" << cmd << "\n";

        std::function<void(const TargetQuery & target_query)> align_read =
            [&hasher, &chainsMutex, &chains_os, &sam_os, &config] (const TargetQuery & target_query) {

                std::optional<chaining::Chain> chain =
                        _map_single(*(std::get<0>(target_query)), *(std::get<1>(target_query)),
                                    config);

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

        std::vector<TargetQuery> targets_queries;
        for (const auto & target : indexed_targets) {
            for (const auto & query : queries) {
                targets_queries.emplace_back(&target, &query);
            }
        }

        process_in_parallel(targets_queries, align_read, nthreads, true);

        chains_os.close();
        sam_os.close();
    }

    void tandem_map(const std::filesystem::path & target_path,
                    const std::filesystem::path & queries_path,
                    const std::filesystem::path & outdir,
                    const bool to_compress,
                    const bool only_index,
                    const size_t nthreads,
                    logging::Logger & logger,
                    const std::string & cmd,
                    const std::filesystem::path & index_path,
                    Config config) {
        if (to_compress) {
            // TODO change that to a parameter call
            StringContig::needs_compressing = true;
        }

        io::SeqReader queries_reader(queries_path);
        std::vector<Contig> queries{queries_reader.readAllContigs()};
        logger.info() << "Queries from " << queries_path << ", total " << queries.size() << " sequences " << std::endl;

        const RollingHash<Config::HashParams::htype> hasher(config.common_params.k, config.hash_params.base);

        io::SeqReader target_reader(target_path);
        std::vector<Contig> targets{target_reader.readAllContigs()};

        for (const Contig & target : targets) {
            logger.info() << "Target length " << target.seq.size() << ", name " << target.id
                          << " from " << target_path << std::endl;
        }

        //const Counter query_counter = kmer_index::_get_readset_counter(hasher, queries);
        //logger.info() << query_counter.size();
        const kmer_index::IndexedContigs indexed_targets =
                kmer_index::get_indexed_targets(queries, targets, outdir, hasher, nthreads, logger, index_path,
                                                config.common_params, config.kmer_indexer_params);

        const auto uncovered_fn = outdir / "norarekmers.bed";
        std::ofstream uncovered_os(uncovered_fn);
        kmer_index::norare_regions2bam(indexed_targets, config.kmer_indexer_params.min_uncovered_len,
                                       config.common_params.k, uncovered_os);
        logger.info() << "Finished exporting long (>= " << config.kmer_indexer_params.min_uncovered_len << " bp) "
                      << "regions without rare k-mers to " << uncovered_fn << std::endl;
        uncovered_os.close();

        if (only_index) {
            if (to_compress) {
                // TODO change that to a parameter call
                StringContig::needs_compressing = false;
            }
            return;
        }

        const auto chains_fn = outdir / "chains.tsv";
        const auto sam_fn = outdir / "alignments.sam";

        logger.info() << "Computing chains and sam records..." << std::endl;
        parallel_run(indexed_targets, queries, hasher, nthreads, chains_fn, sam_fn, cmd, config);

        logger.info() << "Finished outputting chains to " << chains_fn << " and sam records to " << sam_fn << std::endl;

        if (to_compress) {
            // TODO change that to a parameter call
            StringContig::needs_compressing = false;
        }
    }
} // End namespace tandem_mapper