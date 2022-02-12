import subprocess
import sys
from os.path import join

datadir = sys.argv[1]
outdir = sys.argv[2]
veritymap_bin = sys.argv[3]

del_pos=int(sys.argv[4])
del_len=int(sys.argv[5])

additional_option=""
if len(sys.argv) > 6:
    additional_option=" ".join(sys.argv[6:])

DIFF_POS_THRESHOLD=100
DIFF_THRESHOLD=100

fasta_file = join(datadir, "chr9_censat_del2.fasta")
reads_fasta_file = join(datadir, "chr9_censat_reads.fasta")
cmd = veritymap_bin + " --target %s --queries %s -o %s -t 10 %s" % (fasta_file, reads_fasta_file, outdir, additional_option)
print(cmd)
subprocess.call(cmd.split())

maf_file = join(datadir, "chr9_censat_0001.txt")
real_pos = dict()
read_lens=dict()
with open(maf_file) as f:
    for line in f:
        fs = line.split()
        #if "ref" in line:
        #    read_name = "S1_%d" % i
        read_name,ref_s, read_len = fs[0], int(fs[1]), int(fs[2])
        read_lens[read_name] = read_len
        real_pos[read_name] = (int(ref_s),int(read_len))
        #            i += 1

tm_pos = dict()
reads_w_diff = 0
reads_wo_diff = 0
outside_diff = 0

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
        if ref_s >= del_pos:
            shift1 += del_len
        if True: #ref_s+1 < del_pos < ref_e:
            prev_pos, prev_read_pos = 0, 0
            line = f.readline()
            while True:
                line = f.readline()
                if not line or 'Aln' in line: break
                fs = line.split()
                read_pos, ref_pos = int(fs[0]), int(fs[1])
                ref_diff = abs(ref_pos - prev_pos)
                read_diff = abs(read_pos - prev_read_pos)
                diff = abs(ref_diff-read_diff)
                if prev_pos and prev_pos < del_pos < ref_pos:
                    if abs(del_len)-100 < abs(diff) < abs(del_len)+100: 
                        #print(read_name,diff, line)
                        reads_w_diff+=1
                    else:
                        #print("Mapped with incorrect diff", read_name,diff)
                        reads_wo_diff +=1
                elif prev_pos and abs(diff) > DIFF_THRESHOLD:
                    #print("Mapped with diff", diff, "outside the del_len:", read_name, prev_pos)
                    outside_diff +=1
                prev_pos = ref_pos
                prev_read_pos = read_pos
        #read_len = read_lens[read_name]
        read_shift = read_s
        tm_pos[read_name] = (max(0,ref_s-read_shift+shift1),10000)

a=0
b = 0
for read_name in tm_pos:
    if abs(tm_pos[read_name][0] - real_pos[read_name][0]) >= DIFF_POS_THRESHOLD:
        b += 1
        #print("Wrongly mapped", read_name, tm_pos[read_name][0], real_pos[read_name][0])
    else: a+=1
print("Total mapped reads",a+b)
print("Wrongly mapped reads",b)
print("Chains extended through del_len with correct diff", reads_w_diff)
print("Chains extended through del_len with incorrect diff", reads_wo_diff, "(should be 0)")
print("Reads with discrepancies outside the del_len", outside_diff)

MIN_READS_W_DIFF = 5
MAX_READS_WO_DIFF = 0
MAX_OUTSIDE_DIFF = 20

if reads_w_diff <= MIN_READS_W_DIFF or reads_wo_diff > MAX_READS_WO_DIFF or outside_diff > MAX_OUTSIDE_DIFF:
    print(f"Failure: "
          f"MIN_READS_W_DIFF = {MIN_READS_W_DIFF}, real = {reads_w_diff}. "
          f"MAX_READS_WO_DIFF = {MAX_READS_WO_DIFF}, real = {reads_wo_diff}. "
          f"MAX_OUTSIDE_DIFF = {MAX_OUTSIDE_DIFF}, real = {outside_diff}. ")
    sys.exit(1)

print("Successful test on a dataset with difference of %d bp length on %d bp" % (del_len, del_pos))