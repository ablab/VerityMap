from Bio import SeqIO


def get_fasta_lenghts(fasta_fname):
    ref_names = []
    with open(fasta_fname) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
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