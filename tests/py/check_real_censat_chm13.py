import argparse
import os
from pathlib import Path
import subprocess
import sys

from joblib import Parallel, delayed


SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
DEFAULT_INPUT = "/Poppy/abzikadze/centroFlye/centroFlye_repo/data-share/tandemtools2/censat"
DEFAULT_VMBIN = os.path.join(SCRIPT_DIR, os.pardir, os.pardir,
                             'build', 'bin', 'veritymap')
DEFAULT_CENSAT_BED = os.path.join(DEFAULT_INPUT, "cenAnnotation.merged.bed")
DEFAULT_REF_FAI = os.path.join(DEFAULT_INPUT, "ref_v1", "chm13.draft_v1.0.fasta.fai")


class CenSatInfo:
    def __init__(self, chrom, s, e):
        self.chrom = chrom
        self.s = int(s)
        self.e = int(e)

    def __len__(self):
        return self.e - self.e


def get_censats_info(censat_bed_fn, ref_fai_fn):
    censats_info = {}
    with open(censat_bed_fn) as f:
        for line in f:
            chrom, s, e = line.strip().split()
            censats_info[chrom] = CenSatInfo(chrom, s, e)
    with open(ref_fai_fn) as f:
        for line in f:
            chrom, full_len = line.strip().split()[:2]
            full_len = int(full_len)
            if chrom in censats_info:
                censats_info[chrom].full_len = full_len
    return censats_info


def parallel_process(outdir, datadir, censat_info, veritymap_bin, threads,
                     only_index):
    outdir = os.path.join(outdir, censat_info.chrom)
    queries_fn = os.path.join(datadir, "real_reads", "censat",
                              censat_info.chrom + "_censat.fasta")
    target_fn = os.path.join(datadir, "references",
                             censat_info.chrom + "_censat.fasta")
    cmd = veritymap_bin + " --target %s --queries %s -o %s -t %d" % \
        (target_fn, queries_fn, outdir, threads)
    if only_index:
        cmd += " --only-index"
    print(cmd)
    subprocess.call(cmd.split())
    print(censat_info.chrom, " FINISHED!")


def merge_sam(outdir, censats_info, cmd, threads):
    merged_sam_fn = os.path.join(outdir, 'alignments_merged.sam')
    with open(merged_sam_fn, 'w') as f:
        for chrom, censat_info in censats_info.items():
            print(f'@SQ\tSN:{chrom}\tLN:{censat_info.full_len}', file=f)
        print(f'@PG\tID:VerityMap\tPN:VerityMap\tVN:2.0\tCL:{cmd}',
              file=f)
    with open(merged_sam_fn, 'a') as f:
        for chrom, censat_info in censats_info.items():
            print(chrom)
            sam_fn = os.path.join(outdir, censat_info.chrom,
                                  'alignments.sam')
            awk_cmd = \
                f"awk 'NR>2 {{OFS=\"\t\"; $3=\"{chrom}\"; $4+={censat_info.s}; print}}' {sam_fn}"
            subprocess.call(awk_cmd, shell=True, stdout=f)

    merged_bam_fn = os.path.join(outdir, 'alignments_merged.bam')
    sam2bam_cmd = f"samtools view -b {merged_sam_fn} | samtools sort -@{threads} -o {merged_bam_fn}"
    subprocess.call(sam2bam_cmd, shell=True)
    bamindex_cmd = f"samtools index {merged_bam_fn}"
    subprocess.call(bamindex_cmd, shell=True)
    rmsam_cmd = f"rm {merged_sam_fn}"
    subprocess.call(rmsam_cmd, shell=True)


def merge_indexes(outdir, censats_info, cmd):
    merge_indexes_fn = os.path.join(outdir, 'kmer_indexes.tsv')
    outfn = os.path.join(outdir, 'kmer_indexes_merged.tsv')
    with open(outfn, 'w') as fout:
        for chrom, censat_info in censats_info.items():
            fn = os.path.join(outdir, chrom, 'kmer_indexes.tsv')
            with open(fn) as fin:
                for line in fin:
                    _, coord = line.split()
                    coord = int(coord)
                    coord += censats_info[chrom].s
                    fout.write(f'{censats_info[chrom].chrom}\t{coord}\n')


def merge_results(outdir, censats_info, cmd, threads, only_index):
    merge_indexes(outdir, censats_info, cmd)
    if only_index:
        return
    merge_sam(outdir, censats_info, cmd, threads)


def main():
    cmd = ' '.join(sys.argv)
    parser = argparse.ArgumentParser()
    parser.add_argument("--datadir", "-i",
                        default=DEFAULT_INPUT)
    parser.add_argument("--threads", type=int, default=10)
    parser.add_argument("--n-jobs", type=int, default=3)
    parser.add_argument("--outdir", "-o", required=True)
    parser.add_argument("--veritymap-bin",
                        default=DEFAULT_VMBIN)
    parser.add_argument("--censat-bed", default=DEFAULT_CENSAT_BED)
    parser.add_argument("--ref-fai", default=DEFAULT_REF_FAI)
    parser.add_argument("--only-index", action="store_true")
    params = parser.parse_args()

    Path(params.outdir).mkdir(parents=True, exist_ok=True)
    censats_info = get_censats_info(params.censat_bed, params.ref_fai)

    Parallel(n_jobs=params.n_jobs)(delayed(parallel_process)\
        (params.outdir, params.datadir, censat_info,
         params.veritymap_bin, params.threads, params.only_index)
        for censat_info in censats_info.values())

    merge_results(params.outdir, censats_info, cmd,
                  params.threads, params.only_index)


if __name__ == "__main__":
    main()
