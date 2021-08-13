from Bio import SeqIO

def get_fasta_len(fasta_fname):
    with open(fasta_fname) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            return len(record.seq)


def calculate_coverage(assembly_len, bed_fname, read_names=None):
    coverage = [0] * assembly_len
    starts = [0] * assembly_len
    ends = [0] * assembly_len
    with open(bed_fname) as f:
        for line in f:
            fs = line.split()
            ref, ref_s, ref_e, read_name, align_start, align_end, read_len = fs
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
    return coverage