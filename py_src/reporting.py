import base64
from os.path import join
from collections import Counter, defaultdict

import numpy as np

import plotly.graph_objects as go
from plotly.subplots import make_subplots

MIN_COVERAGE = 3
MIN_AF = 20

def format_func(value, tick_number):
    N = value/1000000
    if N == 0:
        return "0"
    else:
        return "%d Mb" % N


def make_plotly_html(assemblies, all_data, out_dir):
    fig = make_subplots(rows=len(all_data), cols=1,
                        subplot_titles=[a.label for a in assemblies])
    step = 200
    for plot_idx, (errors, coverage) in enumerate(all_data):
        asm_id = assemblies[plot_idx].label
        customdata = []
        data = dict()
        data['coverage'] = coverage
        data['coverage'] = [max(data['coverage'][i:i+step]) for i in range(0, len(data['coverage']), step)]
        data['reads'] = [0] * len(data['coverage'])
        data['stddev'] = [0] * len(data['coverage'])
        data['diff'] = [[] for i in range(len(data['coverage']))]
        vals = [0] * len(data['coverage'])
        diffs = [0] * len(data['coverage'])
        reads = [0] * len(data['coverage'])
        for e in errors:
            for i in range(e[0], e[1], step):
                real_pos = i
                real_pos = int(real_pos / step)
                if real_pos >= len(data['reads']):
                    break
                data['reads'][real_pos] += 1
                data['diff'][real_pos].append((e[2], e[3]))
        new_errors =[]
        for i in range(len(data['reads'])):
            diff_arr = [d[1] for d in data['diff'][i]]
            mean_diff = np.median(diff_arr) if diff_arr else 0
            stddev = np.std(diff_arr) if diff_arr else 0
            if stddev > 200:
                stddev = min(stddev, 500)
                filt_reads = [d[0] for d in data['diff'][i] if abs(d[1]-mean_diff) <= min(abs(mean_diff)/2, 3*stddev)]
                filt_diff = [d[1] for d in data['diff'][i] if abs(d[1]-mean_diff) <= min(abs(mean_diff)/2, 3*stddev)]
            else:
                filt_diff = diff_arr
                filt_reads = [d[0] for d in data['diff'][i]]
            reads[i] = filt_reads
            mean_diff2 = np.mean(filt_diff) if filt_diff else 0
            stddev2 = np.std(filt_diff) if filt_diff else 0
            vals[i] = len(filt_reads)*100.0/(data['coverage'][i]+data['reads'][i]) if (data['coverage'][i]+data['reads'][i])>=MIN_COVERAGE and len(filt_reads) > 1 else 0
            diffs[i] = mean_diff2
            customdata.append((len(filt_reads), data['coverage'][i]+data['reads'][i], mean_diff2, stddev2))
            if vals[i] > MIN_AF:
                new_errors.append((i*step, len(filt_reads), data['coverage'][i], data['coverage'][i]+data['reads'][i], mean_diff2, stddev2))

        real_x = [i*step for i in range(len(coverage))]
        fig.add_trace(
            go.Scatter(x=real_x, y=vals, showlegend=False, customdata = customdata,
                       hovertemplate="%{customdata[0]} out of %{customdata[1]} reads, "
                                     "mean diff %{customdata[2]:.2f} std deviation %{customdata[3]:.2f}"), row=plot_idx+1, col=1)
        fig.update_yaxes(range=[-3,105],title_text="% deviated reads", titlefont=dict(size=18), tickfont=dict(size=18),
                         hoverformat="d", row=plot_idx+1, col=1)
        fig.update_xaxes(title_text="Position", titlefont=dict(size=18), tickfont=dict(size=18), hoverformat="d",
                         row=plot_idx+1, col=1)
        bed_fname = join(out_dir, asm_id + "_kmers_dist_diff.bed")
        prev_i = 0
        prev_diff = 0
        with open(bed_fname, "w") as f:
            for x,v,sv_diff in zip(real_x,vals,diffs):
                if v >=20:
                    if not prev_i:
                        prev_i = x
                        prev_v = v
                        prev_diff = sv_diff
                elif prev_i:
                    f.write("%s\t%d\t%d\t%d\t%2.f\n" %
                            (assemblies[plot_idx].contig_name,prev_i, x-step, prev_diff, prev_v))
                    prev_i = 0
                    prev_v = 0
        errors_fname = join(out_dir, asm_id + "_errors.bed")
        with open(errors_fname, "w") as f:
            for e in new_errors:
                f.write(" ".join([str(s) for s in e]))
                f.write("\n")
        reads_fname = join(out_dir, asm_id + "_reads_dist_diff.txt")
        prev_i = 0
        with open(reads_fname, "w") as f:
            support_reads = set()
            for x,v,r in zip(real_x,vals,reads):
                if v >=20:
                    if not prev_i:
                        prev_i = x
                        prev_v = v
                    for read in r:
                        support_reads.add(read)
                elif prev_i:
                    f.write("%s\t%d\t%d\n" %
                            (assemblies[plot_idx].contig_name,prev_i, x-step))
                    for r in support_reads:
                        f.write(r + "\n")
                    prev_i = 0
                    support_reads = set()
        plot_idx += 1
    plot_fname = join(out_dir, "kmers_dist_diff.html")
    fig.write_html(plot_fname)
    print("  Difference in k-mer distances plot saved to %s" % plot_fname)

