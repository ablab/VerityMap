from Bio import SeqIO
import gzip
import os


def get_asm_lenghts(fname):
    ref_names = []

    opener, fname_wogzip = (gzip.open, os.path.splitext(fname)[0]) \
        if fname.endswith('.gz') else (open, fname)
    _, ext = os.path.splitext(fname_wogzip)
    ext = ext[1:]
    if ext in ['fa', 'fasta']:
        formt = 'fasta'
    elif ext in ['fq', 'fastq']:
        formt = 'fastq'
    else:
        raise ValueError("Can't guess format of " + fname +
                         " from its extension " + ext)

    with opener(fname, 'rt') as handle:
        for record in SeqIO.parse(handle, formt):
            ref_names.append((record.id, len(record.seq)))
    return ref_names


def calculate_coverage(assembly_lenghts, bed_fname, read_names=None):
    asm_cov = dict()
    for ref_name, assembly_len in assembly_lenghts:
        coverage = [0] * assembly_len
        starts = [0] * assembly_len
        ends = [0] * assembly_len
        with open(bed_fname) as f:
            for line in f:
                fs = line.split()
                ref, ref_s, ref_e, read_name, align_start, align_end, read_len = fs
                if ref != ref_name: continue
                if read_names is not None and read_name not in read_names:
                    continue
                ref_s, ref_e, align_start, align_end, read_len = map(int, (ref_s, ref_e, align_start, align_end, read_len))
                starts[max(0,ref_s)] += 1
                ends[min(assembly_len-1, ref_e - 1)] += 1
        cur_cov = 0
        for i in range(assembly_len):
            cur_cov += starts[i]
            cur_cov -= ends[i]
            coverage[i] = cur_cov
        asm_cov[ref_name] = coverage
    return asm_cov
