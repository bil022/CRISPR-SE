# CRISPR-SE
CRISPR-SE uses fasta format: (ex: mm10.fa, hg38.fa, ecoli.fa, simple.fa). There are two options: 1. standard fasta files 2. simple format, list of 20-bp guide RNAs (gRNAs), one gRNA per line, for example:
```
>simple.fa
TCTATTTTGTGGTTACTTTG
GTGGTTACTTTGAGGAGAGT
CTAAATCAGGATCAGATTCA
```
## Compile CRISPR-SE
```
git clone https://github.com/bil022/CRISPR-SE
cd CRISPR-SE
make
./se
```
## Usage
```
Crispr-SE <command> [options]

Command:
  --index  create index from reference/query sequence in the FASTA format
  --build  build whole genome single guide RNA (gRNA)

Options:
  -p INT  number of threads [2]
  -r STR  reference genome id (mm9, mm10, hg19, hg38, etc)
  -s  The FASTA format is simple format of 20-nt gRNA per line
  -q  The query (user inputs, search reference genome if not set)
  -m INT  Max mismatch, 0 for CREST-SE(2*#seed+#distal<4), 1 or more for #mismatches, default: 0
  -n INT  Max off-target, 0 for all, default: 1
  -v  verbose mode, default: off
```
## Create index
#### Create index for ecoli reference genome
```
./se --index -r ecoli
```
#### Create index for simple.fa with simple format (one gRNA per line)
```
./se --index -r simple -s
```
## Search gRNAs
#### Search genome-wide gRNA for ecoli
```
./se --build -r ecoli -p 2
```
#### Search gRNAs in simple.fa using ecoli as reference genome
```
./se --build -r ecoli -q simple -p 2
```
#### Dump off-targets for gRNAs in simple.fa
```
./se --build -r ecoli -q simple -p 2 -v

```
#### Convert into bam format
For reference genome only
```
cat ecoli.h ecoli.mm | samtools view -Sb - > ecoli.bam
```
For user input simple.fa
```
cat simple.h simple.mm | samtools view -Sb - > simple.bam
```
#### View results
For reference genome
```
samtools view ecoli.bam | less
```
For user input simple.fa
```
samtools view ecoli.bam | less
```
