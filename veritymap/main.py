import datetime
import os
import sys
from os.path import isdir, abspath

root = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)
sys.path.insert(0, root)

import click as click

from veritymap.py_src.assembly import Assembly
from veritymap.py_src.mapper import do


@click.command()
@click.argument('assembly_fnames', type=click.Path(exists=True), nargs=-1)
@click.option('--reads', 'reads_fname', type=click.Path(exists=True), help='File with ONT/PacBio reads')
@click.option('-o', 'out_dir', type=click.Path(), required=True, help='Output folder')
@click.option('-t', 'threads', type=click.INT, help='Threads', default=4)
@click.option('-d', 'datatype', type=click.Choice(['hifi', 'ont']), help='Sequencing platform, supported types are: '
                                                                         '"hifi" for PacBio HiFi reads and "ont" for ONT reads.'
                                                                         'Please note that "ont" mode is experimental and '
                                                                         'should be used with extra care')
@click.option('-f', '--no-reuse', 'no_reuse', is_flag=True, help='Do not reuse old files')
@click.option('--careful', 'is_careful', is_flag=True, help='Run mapper in a careful mode to better detect inconsistencies. Can be time- and memory-consuming. Not recommended to run on the whole genome. ')
@click.option('-l', 'labels', help='Comma separated list of assembly labels')
def main(assembly_fnames, reads_fname, labels, out_dir, threads, no_reuse, is_careful, datatype):
    date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("%s VerityMap started" % date)
    if not reads_fname:
        print("ERROR! You should specify ONE path to a file with reads (ONT or Pacbio CLR reads)")
        sys.exit(2)

    if not assembly_fnames:
        print("ERROR! You should specify at least one assembly file.")
        sys.exit(2)

    if not datatype:
        datatype = "hifi"

    list_labels = [None] * len(assembly_fnames)
    if labels:
        list_labels = labels.replace('"', '').split(',')
        if len(list_labels) != len(assembly_fnames):
            print("ERROR! Number of labels must correspond to the number of analyzed assemblies")
            sys.exit(2)

    assemblies = [Assembly(assembly_fnames[i], name=list_labels[i], out_dir=out_dir) for i in range(len(assembly_fnames))]
    out_dir = abspath(out_dir)
    if not isdir(out_dir):
        os.makedirs(out_dir)
    do(assemblies, reads_fname, datatype, out_dir, threads, no_reuse, is_careful)
    date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print()
    print("%s VerityMap finished!" % date)


if __name__ == '__main__':
    main()
