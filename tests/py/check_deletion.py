import subprocess
import sys
from os.path import join

point_del = 456400
deletion = 10000
threshold=100

real_pos = dict()
i = 1
read_lens=dict()

datadir = sys.argv[1]
outdir = sys.argv[2]
veritymap_bin = sys.argv[3]

fasta_file = join(datadir, "chr7_ext_del.fasta")
reads_fasta_file = join(datadir, "chr7_ext_subreads.fasta")
cmd = veritymap_bin + " --target %s --queries %s -o %s -t 10" % (fasta_file, reads_fasta_file, outdir)
print(cmd)
subprocess.call(cmd.split())

maf_file = join(datadir, "chr7_ext_0001.maf")
with open(maf_file) as f:
    for line in f:
        fs = line.split()
        if "ref" in line:
            read_name = "S1_%d" % i
            read_len,ref_s = int(fs[3]), int(fs[2])
            read_lens[read_name] = read_len
            real_pos[read_name] = (int(ref_s),int(read_len))
            i += 1

tm_pos = dict()
reads_w_diff = 0
reads_wo_diff = 0

chains = join(outdir, "chains.tsv")
with open(chains) as f:
    for line in f:
        fs = line.split()
        if 'Aln' not in fs[0]: continue
        read_name, ref_name, read_s, read_e, read_len, ref_s, ref_e = fs[1:8]
        read_name=read_name.replace('+','').replace('-','')
        try: ref_s, ref_e, read_s, read_e = map(int, (ref_s, ref_e, read_s, read_e))
        except: continue
        shift1 = 0
        if ref_s >= point_del:
            shift1 += deletion
        if ref_s < point_del and ref_e > point_del:
            prev_pos, prev_read_pos = 0, 0
            line = f.readline()
            while True:
                line = f.readline()
                if not line or 'Aln' in line: break
                fs = line.split()
                read_pos, ref_pos = int(fs[0]), int(fs[1])
                if prev_pos <= point_del and ref_pos > point_del:
                    ref_diff = abs(ref_pos - prev_pos)
                    read_diff = abs(read_pos - prev_read_pos)
                    diff = ref_diff-read_diff
                    if -10500 < diff < -9500: reads_w_diff+=1
                    else: reads_wo_diff +=1
                prev_pos = ref_pos
                prev_read_pos = read_pos
        read_len = read_lens[read_name]
        read_shift = read_s
        tm_pos[read_name] = (max(0,ref_s-read_shift+shift1),10000)

a=0
b = 0
for read_name in tm_pos:
    if abs(tm_pos[read_name][0] - real_pos[read_name][0]) >= threshold:
        b += 1
        print("Wrongly mapped", read_name, tm_pos[read_name][0], real_pos[read_name][0])
    else: a+=1
print("Wrongly mapped reads",b, "total mapped", a+b)
print("Chains extended through deletion with 10kbp diff", reads_w_diff, "(prev result: 8)")
print("Chains extended through deletion with incorrect diff", reads_wo_diff, "(should be 0)")

MIN_READS_W_DIFF = 5
MAX_READS_WO_DIFF = 0

if reads_w_diff >= MIN_READS_W_DIFF and reads_wo_diff > MAX_READS_WO_DIFF:
    print(f"Failure: "
          "MIN_READS_W_DIFF = {MIN_READS_W_DIFF}, real = {reads_w_diff}. "
          "MAX_READS_WO_DIFF = {MAX_READS_WO_DIFF}, real = {reads_wo_diff}")
    sys.exit(1)

print("Successful test on a dataset with deletion")
