import argparse
from collections import defaultdict
import re
from itertools import groupby


def get_ref_len(sam_fn):
    ref_len = {}
    with open(sam_fn) as f:
        for line in f:
            if line[0] != '@':
                break
            search = re.search('^@SQ\tSN:(.+)\tLN:(.+)', line)
            if search is None:
                continue
            ref, length = search.group(1), int(search.group(2))
            ref_len[ref] = length
    return ref_len


class Cigar:
    op_chars = ("S", "M", "I", "D")

    def __init__(self, line):
        self.ops = []

        cig_iter = groupby(line, lambda chr: chr.isdigit())
        for _, length_digits in cig_iter:
            length = int(''.join(length_digits))
            op = next(next(cig_iter)[1])
            self.ops.append((op, length))

    def soft_clip_left(self):
        if len(self.ops) and self.ops[0][0] == "S":
            return self.ops[0][1]
        return 0

    def soft_clip_right(self):
        if len(self.ops) and self.ops[-1][0] == "S":
            return self.ops[-1][1]
        return 0

    def ref_len(self):
        length = 0
        for op, op_len in self.ops:
            assert op in self.op_chars
            if op in "MD":
                length += op_len
        return length

    def query_len(self):
        length = 0
        for op, op_len in self.ops:
            if op in "SMI":
                length += op_len
        return length


def get_coverage_nclipped(ref_len, sam_fn, min_clip=1000):
    coverage = {ref: [0] * (length+1) for ref, length in ref_len.items()}
    n_clipped = {ref: defaultdict(int) for ref in ref_len}

    with open(sam_fn) as f:
        while f:
            line = f.readline()
            if line[0] != '@':
                break
        while f and len(line):
            line = line.strip().split('\t')
            if len(line) <= 3:
                print(line)
            ref = line[2]
            ref_s = int(line[3])
            cigar = Cigar(line[5])
            ref_e = ref_s + cigar.ref_len()

            if cigar.soft_clip_left() > min_clip:
                n_clipped[ref][ref_s] += 1
            if cigar.soft_clip_right() > min_clip:
                n_clipped[ref][ref_e-1] += 1

            coverage[ref][ref_s] += 1
            coverage[ref][ref_e] -= 1
            if line[0] == "S2_1117":
                print(ref_s, ref_e, cigar.soft_clip_right())

            line = f.readline()

    for ref, cov_ref in coverage.items():
        for i in range(1, len(cov_ref)):
            cov_ref[i] += cov_ref[i-1]

    return coverage, n_clipped


def get_substantial_clip(coverage, n_clipped, min_ratio=0.3, min_clip=3):
    subst_clipped = {ref: defaultdict(int) for ref in coverage}
    for ref, n_clipped_ref in n_clipped.items():
        cov_ref = coverage[ref]
        for pos, cnt in n_clipped_ref.items():
            if cov_ref[pos] and \
                    cnt / cov_ref[pos] > min_ratio and \
                    cnt > min_clip:
                subst_clipped[ref][pos] += cnt
    return subst_clipped


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--sam", required=True, help="Input sam file")

    params = parser.parse_args()
    ref_len = get_ref_len(params.sam)
    coverage, n_clipped = get_coverage_nclipped(ref_len, params.sam)
    subst_clip = get_substantial_clip(coverage, n_clipped)

    print(subst_clip)


if __name__ == "__main__":
    main()
