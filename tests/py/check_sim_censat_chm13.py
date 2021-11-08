import os
import sys
import subprocess
import pysam
from joblib import Parallel, delayed

parallel_args = {'n_jobs': 3}
threshold=100
threads = 10

datadir=sys.argv[1]
outdir=sys.argv[2]

from pathlib import Path
Path(outdir).mkdir(parents=True, exist_ok=True)

veritymap_bin = sys.argv[3]

chroms = []
chrom_len= dict()
with open(os.path.join(datadir, "cenAnnotation.merged.bed")) as f:
    for line in f:
        chrom, s, e = line.split()
        chr_s, chr_e = int(s), int(e)
        chroms.append((chrom, chr_s, chr_e))
        chrom_len[chrom] = chr_e - chr_s


def parallel_process(outdir, chrom, s, e):
    outdir = outdir + "/" + chrom
    reads_fasta_file = datadir + "/simulated_reads/" + chrom + "_censat_0001.fasta"
    fasta_file = datadir + "/references/" + chrom + "_censat.fasta"
    cmd = veritymap_bin + " --target %s --queries %s -o %s -t %d" % (fasta_file, reads_fasta_file, outdir, threads)
    print(cmd)
    subprocess.call(cmd.split())
    print(chrom, " FINISHED!")

results = Parallel(**parallel_args)(delayed(parallel_process)(outdir, chrom, s, e) for chrom, s, e in chroms)

outfn = outdir + "/sim_results.tsv"
outf = open(outfn, "w")
outf.write("total\tmapped\tmapped_freq\twrong\twrong_perc\tunmapped\tunmapped_perc\tunmapped_bases\tunmapped_bases_perc\tassembly_length\tuncovered_bases\tuncovered_bases_perc\tno_rare_kmers_bases")

def group(number):
    s = '%d' % number
    groups = []
    while s and s[-1].isdigit():
        groups.append(s[-3:])
        s = s[:-3]
    return s + ','.join(reversed(groups))

for c in list(range(1,23))+['X']:
#for c in ['X']:
    real_pos = dict()
    read_lens=dict()
    chr_id = "chr" + str(c)
    i = 1
    total_read_len = 0
    mapped_read_len = 0
    if not os.path.exists("%s/%s/alignments.sam" % (outdir,chr_id)): continue
    with open("%s/simulated_reads/%s_censat_0001.txt" % (datadir, chr_id)) as f:
        for line in f:
            read_name, ref_s, read_len = line.split()
            real_pos[read_name] = [int(ref_s)]
            total_read_len += int(read_len)

    tm_pos = dict()
    assembly_len = chrom_len[chr_id]
    coverage = [0] * assembly_len
    starts = [0] * assembly_len
    ends = [0] * assembly_len
    samfile = "%s/%s/alignments.sam" % (outdir,chr_id)
    for r in pysam.AlignmentFile(samfile, "r").fetch():
        tm_pos[r.query_name] = max(0, r.reference_start - r.query_alignment_start)
        mapped_read_len += r.query_alignment_length
        starts[r.reference_start] += 1
        ends[r.reference_end] += 1

    norarekmers_len = 0
    with open("%s/%s/norarekmers.bed" % (outdir,chr_id)) as f:
        for line in f:
            fs = line.split()
            norarekmers_len += int(fs[3])
    a=0
    b = 0
    for read_name in tm_pos:
        if read_name not in real_pos: continue
        if abs(tm_pos[read_name] - real_pos[read_name][0]) >= threshold:
            b += 1
            # write erroneously mapped reads
            #print(read_name, tm_pos[read_name], real_pos[read_name][0])
        else: a+=1
    cur_cov = 0
    for i in range(assembly_len):
        cur_cov += starts[i]
        cur_cov -= ends[i]
        coverage[i] = cur_cov
    uncovered = sum([1 for c in coverage if c == 0])
    b_pct = (b*100/a)
    total = len(real_pos)
    unmapped = len(real_pos) - a
    unm_pct=unmapped*100/total
    unmapped_read_len = total_read_len - mapped_read_len
    unm_len_pct=unmapped_read_len*100/total_read_len
    num_n = {'chr15': 2700000, 'chr14': 675000, 'chr13': 4050000, 'chr21': 3150000, 'chr22': 900000}
    real_len = 0
    real_uncovered = 0
    if chr_id in num_n:
        assembly_len -= num_n[chr_id]
        uncovered -= num_n[chr_id]
    print("%s. Total: %d, mapped: %d (%.2f%%), wrong: %d (%.2f%%), unmapped: %d reads (%.2f%%) %d bp (%.2f%%). Assembly length: %d bp %s, uncovered bases: %d bp (%.2f%%) %s, regions with no rare k-mers %d bp" %
      (chr_id, total, a, (a*100/total), b, (b*100/total), unmapped, unm_pct, unmapped_read_len, unm_len_pct, assembly_len, ("" if not real_len else "without N's: %d" % real_len), \
       uncovered, uncovered*100/assembly_len, ("" if not real_len else "without N's: %d (%.2f%%)" % (real_uncovered, real_uncovered*100/real_len)), norarekmers_len))
    outf.write("\n%s \t%s \t%s \t%.2f \t%s \t%.2f \t%s \t%.2f \t%s \t%.2f \t%s \t%s \t%.2f \t%d" %
      (chr_id, total, a, (a*100/total), b, (b*100/total), unmapped, unm_pct, unmapped_read_len, unm_len_pct, assembly_len, \
       uncovered, uncovered*100/assembly_len, norarekmers_len))

outf.close()

fail = False
if len(sys.argv) > 4:
    import pandas as pd

    results_df = pd.read_csv(outfn, sep='\t', index_col=0)
    worst = pd.read_csv(sys.argv[4], sep='\t', index_col=0)

    low_mapped_freq = results_df['mapped_freq'][results_df['mapped_freq'] < worst['mapped_freq']]
    if len(low_mapped_freq):
        print('ERROR: low mapped freq:')
        print(low_mapped_freq)
        fail = True
    wrong = results_df['wrong'][results_df['wrong'] > worst['wrong']]
    if len(wrong):
        print('ERROR: too many wrong mappings:')
        print(wrong)
        fail = True
    unmapped_perc = results_df['unmapped_perc'][results_df['unmapped_perc'] > worst['unmapped_perc']]
    if len(unmapped_perc):
        print('ERROR: too many unmapped reads:')
        print(unmapped_perc)
        fail = True
    uncovered_bases_perc = results_df['uncovered_bases_perc'][results_df['uncovered_bases_perc'] > worst['uncovered_bases_perc']]
    if len(uncovered_bases_perc):
        print('ERROR: too many uncovered bases:')
        print(uncovered_bases_perc)
        fail = True

if fail:
    sys.exit(1)

print("Successful run on simulated reads from censat regions")
