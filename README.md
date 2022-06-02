# VerityMap

## Quick start
```shell
    veritymap --reads test_dataset/test_query.fasta test_dataset/test_target.fasta -o test_outdir -d hifi
```

## Introduction

**VerityMap** (formerly known as TandemMapper2) is designed for mapping long reads (PacBio HiFi or ONT) to assemblies of extra-long tandem repeats, such as centromeres, but can be applied to whole-genome assemblies. The tool outputs SAM file that can be used in any downstream analysis. In addition, VerityMap yields an information about possible errors and heterozygous sites in the assembly based on analysis of rare k-mers.

## Installation

Requirements are listed in ```requirements.txt``` and can be installed through Conda as ```conda install --file requirements.txt``` or pip as ```pip install -r requirements.txt```.

## Usage

```shell
veritymap [options] --reads <reads_file> -d <hifi,ont> -o <output_dir> <assembly_file1> <assembly_file2>

Required arguments:
    --reads         PATH                 File with Oxford Nanopore or PacBio HiFi reads used for ETR assembly
    -o              PATH                 Folder to store all result files
    -d              <hifi,ont>           Type of used sequencing platform ("hifi" for PacBio HiFi reads and "ont" for ONT reads)

Optional arguments:
    -t             INT                  Maximum number of threads [default: 4]
    -l             \"label,label,...\"  Human-readable names of assemblies to use in reports, comma-separated. If contain spaces, use quotes
```
In case VerityMap is built locally, it should be run as follows:

```shell
python veritymap/main.py [options] --reads <reads_file> -d <hifi,ont> -o <output_dir> <assembly_file1> <assembly_file2>
```

## Output files

The following files are contained in `<output_dir>` directory (specified by `-o`) and include results
for all input assemblies. 

`<output_dir>/*_alignment.bed` - VerityMap alignments in BED format

`<output_dir>/*_alignment.sam` - VerityMap alignments in SAM format

`<output_dir>/*_kmers_dist_diff.bed` - BED file with coordinates of possible heterozygous sites and errors. The format is  `{ref_name} {start} {end} {misassembly_len} {% discordant reads}`. Please do not pay much attention to the last value, and use the interactive HTML plot described below to see % discordant reads. Negative misassembly length corresponds to the deletion in the reference.
`<output_dir>/*_kmers_dist_diff.html` - interactive HTML plot showing possible heterozygous sites and errors.
VerityMap analyzes distances between consecutive rare k-mers in the assembly and in a read. If these distances are inconsistent, i.e., a read does not support the assembly in this position, it may indicate the presence of a heterozygous site (if the percent of deviated reads is 20-80%) or an assembly error (if the percent of such reads is 80-100%). 
In the plot, OX shows position in the assembly, OY shows the percent of reads that disagree with the assembly in this position. By hovering over the plot, you can also get the information about the number of supporting reads, the mean difference in distances (essentially, the length of an indel) and standard deviation of the difference.
In the BED file, VerityMap outputs approximate coordinates, length, and frequency of detected variants. Importantly, reported coordinates are coordinates of the rare k-mers that are nearest to the variant, i.e., the real coordinates of the variant are unknown and locate somewhere inside the reported region. You can find the real coordinates analyzing the SAM file with the alignments. 

The remaining output files are technical and likely should not concern the end user.

## Citation

Currently, you can refer to the first TandemTools paper:

Alla Mikheenko, Andrey V Bzikadze, Alexey Gurevich, Karen H Miga, Pavel A Pevzner, TandemTools: mapping long reads and assessing/improving assembly quality in extra-long tandem repeats, Bioinformatics, Volume 36, Issue Supplement_1, July 2020, Pages i75–i83, https://doi.org/10.1093/bioinformatics/btaa440

The paper describing VerityMap algorithm is in preparation.

## Contacts

Please report any problems to the [issue tracker](https://github.com/ablab/tandemQUAST/issues). Alternatively, you can write directly to [a.mikheenko@spbu.ru](mailto:a.mikheenko@spbu.ru).
