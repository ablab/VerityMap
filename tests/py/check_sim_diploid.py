import os
import sys


threshold = 100
chains_file = sys.argv[1]
ANSWER_DIR = sys.argv[2]

real_pos = dict()
i = 1
read_scores = dict()
with open(os.path.join(ANSWER_DIR, "maf_first.txt")) as f:
    for line in f:
        read_name = "S1_%d" % i
        read_name, ref_s = line.split()[:2]
        real_pos[read_name] = int(ref_s)
        i += 1
i = 1
with open(os.path.join(ANSWER_DIR, "maf_second.txt")) as f:
    for line in f:
        read_name, ref_s = line.split()[:2]
        read_name = "S2_%d" % i
        real_pos[read_name] = int(ref_s)
        i += 1

tm_pos = dict()
d = 0
with open(chains_file) as f:
    for line in f:
        fs = line.split()
        if 'Aln' not in fs[0]:
            continue
        read_name, ref_name, read_s, read_e, read_len, ref_s, ref_e = fs[1:8]
        if 'S2' in read_name and 'hg' not in ref_name:
            d += 1
            print(fs)
            continue
        if 'S2' not in read_name and 'hg' in ref_name:
            d += 1
            print(fs)
            continue
        score = float(fs[-1])
        read_name = read_name.replace('+', '').replace('-', '')
        try:
            ref_s, ref_e, read_s, read_e = map(
                int, (ref_s, ref_e, read_s, read_e))
        except BaseException:
            continue
        shift1 = 0
        read_shift = read_s  # if strand != 'backward' else read_len - read_e
        if read_name in read_scores and read_scores[read_name] > score:
            continue
        read_scores[read_name] = score
        tm_pos[read_name] = (max(0, ref_s-read_shift+shift1), 10000)

a = 0
b = 0
for read_name in tm_pos:
    if abs(tm_pos[read_name][0] - real_pos[read_name]) >= threshold:
        b += 1
        #if b <10:print(read_name, tm_pos[read_name][0], real_pos[read_name])
    else:
        #if a <4:print(read_name, tm_pos[read_name][0], real_pos[read_name])
        a += 1
print("Wrong target: %d reads" % d)
print("Wrong mapped:", b, "reads, correct mapped:",
      a, "reads, unmapped:", len(real_pos) - (a+b))


if len(sys.argv) >= 4:
    outfn = sys.argv[3]
    with open(outfn, 'w') as f:
        print("Wrong_target\t%d" % d, file=f)
        print("Misplaced\t", b, "\nCorrect\t", a, "\nUnmapped\t",
              len(real_pos) - (a+b), sep='', file=f)
