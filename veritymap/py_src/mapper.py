import os
import subprocess
import sys
import shutil
from collections import defaultdict
from os.path import abspath, exists, isdir, join

from veritymap.py_src.reporting import *
from veritymap.py_src.utils import *

source_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
MAPPER_BIN = join(source_dir, "build/bin/veritymap")
MIN_CHAIN_KMERS = 5
MIN_CHAIN_LEN = 2000
MAX_DIFF_PB = 100
MAX_DIFF_ONT = 1000


def make_mapper():
    if exists(MAPPER_BIN):
        return

    print('Compiling mapper...')
    return_code = subprocess.call(['make', '-C', source_dir], stdout=open("make.log","w"), stderr=open("make.err","w"))
    if return_code != 0:
        raise
    print('Mapper is compiled successful!')


def run_mapper(assembly, reads_fname, out_dir, threads, datatype, is_careful):
    try:
        make_mapper()
    except:
        print('Failed to compile VerityMap! Please try to compile it manually: run "make" in %s'
              % source_dir)
        sys.exit(2)

    cmd = [MAPPER_BIN,
           '--target', assembly.fname, '--queries', reads_fname,
           '-o', join(out_dir, 'veritymap'), '-t', str(threads), '--config', datatype]
    if is_careful: cmd += ['--careful']
    subprocess.call(cmd)
    output_fname = join(out_dir, 'veritymap', 'chains.tsv')
    sam_fname = join(out_dir, 'veritymap', 'alignments.sam')
    shutil.move(output_fname, assembly.chains_fname)
    shutil.move(sam_fname, assembly.sam_fname)
    print('Mapping finished!')


def postprocess_chains(assembly, datatype):
    read_alignments = defaultdict(list)
    read_lengths = dict()
    read_seeds = defaultdict(lambda: defaultdict(list))
    with open(assembly.chains_fname) as f:
        for line in f:
            fs = line.split()
            if "Aln" in line and len(fs) >= 8:
                read_name, ref_name, align_start, align_end, read_len, ref_start, ref_end = fs[1:8]
                align_start, align_end, read_len, ref_start, ref_end = map(int, (align_start, align_end, read_len, ref_start, ref_end))
                read_name = read_name.replace('-','')
                read_lengths[read_name] = read_len
                read_alignments[read_name].append((ref_name, ref_start, ref_end, align_start, align_end))
            elif "Aln" not in line and len(fs) >= 2:
                read_pos, ref_pos = int(fs[0]), int(fs[1])
                read_seeds[read_name][(ref_name, ref_start, ref_end, align_start, align_end)].append((read_pos, ref_pos))
    num_alignments = 0
    all_errors = defaultdict(list)
    max_diff = MAX_DIFF_ONT if datatype == 'ont' else MAX_DIFF_PB
    with open(assembly.bed_fname, "w") as f:
        for read_name, aligns in read_alignments.items():
            max_kmers = 0
            selected_chain = []
            selected_errors = []

            for a_i,align in enumerate(aligns):
                seeds = read_seeds[read_name][align]
                seeds.sort(key=lambda x: x[1])
                best_chain = None
                best_kmers = 0
                best_errors = []
                cur_errors = []
                ref_name = align[0]
                if len(seeds) >= MIN_CHAIN_KMERS:
                    new_chains = []
                    breakpoints = []
                    for i in range(1, len(seeds)):
                        ref_s, ref_e, aln_s, aln_e = seeds[i-1][1],seeds[i][1],seeds[i-1][0],seeds[i][0]
                        ref_diff = abs(ref_e - ref_s)
                        read_diff = abs(aln_e - aln_s)
                        if abs(ref_diff-read_diff) >= max_diff:
                            breakpoints.append(i-1)
                            cur_errors.append((seeds[i-1][1], seeds[i][1], read_name, ref_diff-read_diff))
                    if breakpoints:
                        chain_start1, chain_end1, chain_start2, chain_end2 = seeds[0][1], seeds[-1][1], seeds[0][0], seeds[-1][0]
                        start_n = 0
                        for p in breakpoints:
                            chain_end1, chain_end2 = seeds[p][1], seeds[p][0]
                            new_chains.append((chain_start1, chain_end1, chain_start2, chain_end2, p-start_n+1))
                            if p < len(seeds):
                                chain_start1, chain_start2, start_n = seeds[p+1][1], seeds[p+1][0], p+1

                        chain_end1, chain_end2 = seeds[-1][1], seeds[-1][0]
                        if chain_end1 > chain_start1:
                            new_chains.append((chain_start1, chain_end1, chain_start2, chain_end2, len(seeds) - p))
                        chains = []
                        total_kmers = 0
                        total_len = 0
                        for c in new_chains:
                            chain_start1, chain_end1, chain_start2, chain_end2, chain_kmers = c
                            if chain_kmers > MIN_CHAIN_KMERS and chain_end1 - chain_start1 >= MIN_CHAIN_LEN:
                                chains.append([chain_start1, chain_end1, chain_start2, chain_end2])
                                total_kmers += chain_kmers
                                total_len += chain_end1-chain_start1
                        if total_kmers > best_kmers:
                            best_kmers = total_kmers
                            best_chain = chains
                            best_errors = cur_errors
                    else:
                        best_kmers = len(seeds) if len(seeds) > MIN_CHAIN_KMERS/2 else 0
                        best_chain = [[seeds[0][1], seeds[-1][1], seeds[0][0], seeds[-1][0]]]
                        best_errors = []

                    if best_kmers > max_kmers:
                        max_kmers = best_kmers
                        selected_chain = best_chain
                        selected_errors = best_errors
                        best_ref = ref_name

            for c in selected_chain:
                ref_start, ref_end, align_start, align_end = c
                if (ref_end-ref_start) < MIN_CHAIN_LEN:
                    continue
                num_alignments += 1
                f.write("%s\t%d\t%d\t%s\t%d\t%d\t%d\n" %
                        (best_ref, ref_start, ref_end, read_name, align_start, align_end, read_lengths[read_name]))
            all_errors[ref_name].extend(selected_errors)

    print("  Total %d alignments" % num_alignments)
    print("  Longest chains saved to %s" % assembly.bed_fname)
    print("  SAM file with alignments saved to %s" % assembly.sam_fname)
    return all_errors


def do(assemblies, reads_fname, datatype, out_dir, threads, no_reuse, is_careful):
    print("")
    print("*********************************")
    print("Read mapping started...")
    assemblies_to_process = [assembly for assembly in assemblies if not exists(assembly.bed_fname) or no_reuse]
    for assembly in assemblies_to_process:
        run_mapper(assembly, reads_fname, out_dir, threads, datatype, is_careful)
    all_data = []
    for assembly in assemblies:
        errors = postprocess_chains(assembly, datatype)
        coverage = calculate_coverage(get_fasta_lenghts(assembly.fname), assembly.bed_fname)
        all_data.append((errors, coverage))
    make_plotly_html(assemblies, all_data, out_dir)
