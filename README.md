# CRISPR-SE
## Input:
Genome files in fasta format: (ex: mm10.fa, hg38.fa, ecoli.fa, simple.fa). There are two options: 1. standard fasta files 2. simple format, list of 20-bp guide RNAs (gRNAs), one gRNA per line, for example:
```bash
>simple
TCTATTTTGTGGTTACTTTG
GTGGTTACTTTGAGGAGAGT
CTAAATCAGGATCAGATTCA
```
## Compile CRISPR-SE
```
$ git clone https://github.com/bil022/CRISPR-SE
$ cd CRISPR-SE
$ make
$ ./se
Program: Crispr-SE (CRISPR Search Engine)
Contact:  Bin Li <bil022@ucsd.edu>
Usage:  Crispr-SE <command> [options]

Command:
--index  create index from reference/query sequence in the FASTA format
--build  build whole genome single guide RNA (gRNA)

Options:
-p INT  number of threads [2]
-r STR  reference genome id (mm9, mm10, hg19, hg38, etc)
-s  The FASTA format is simple format of 20-nt gRNA
-q  The indexed query
-m INT  Max mismatch, 0 for CREST-Seq, 1+ for #mismatches, default: 0
-n INT  Max off-target, 0 for all, default: 1
-v  verbose mode, default: false
```
## Create index
```
# create index for ecoli reference genome
$ ./se --index -r ecoli
# create index for simple.fa with one gRNA per line
$ ./se --index -sr simple
```
## Search gRNAs
```
# search genome-wide gRNA
$ ./se --build -r ecoli -p 2

# search gRNA in list of gRNAs in simple format
$ ./se --build -r ecoli -q simple -p 2

# search offtargets (slower)
$ ./se --build -r ecoli -q simple -p 2 -v

```
