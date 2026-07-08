# Landrace_v1
Scripts used for the creation of the Landrace_v1 assembly.

## Assembly

### Flye

Reads were filtered using Filtlong
```
filtlong --min_length 4000 --keep_percent 90 Landrace_V1.ONT_rawreads.fq.gz | pigz -p4 > Landrace_V1.ONT_rawreads.filtered.fq.gz
```

Filtered ONT reads were assembled into contigs using Flye (v 2.9.1-b1780) 
```
flye \
  --nano-raw Landrace_V1.ONT_rawreads.filtered.fq.gz \
  --min-overlap 20000 \
  --genome-size 2.5g \
  --out-dir Landrace_V1_20k_flye \
  --asm-coverage 40 \
  -t 30
```
The `--asm-coverage` flag was set to 40 to reduce the memory consumption in the initial disjoining step
A range of required read overlaps (`--min-overlap` flag) were used (5k, 7k, 10k, 15k, 20k, 25k and 30k) to produce multiple independent assemblies. 

### HiFiasm

Raw reads from PacBio HiFi were filtered with hifiadapterfilt
```
bash hifiadapterfilt.sh -p Landrace_V1.hifi_rawreads -t 30
```

Main assembly was created with filtered ONT reads and filtered HiFi reads using HiFiasm (v 0.19.6-r595) with default settings
```
hifiasm \
  -o Landrace_hifiasm \
  -t 30 \
  --ul Landrace_V1.ONT_rawreads.filtered.fq.gz \
  Landrace_V1.hifi_rawreads.filt.fastq.gz
```

### Merging
Align HiFi assembly contigs against Flye assembly contigs.
```
<software> Landrace_V1_20k.flye.fasta Landrace_V1.hifi.fasta
```

Identify possible cases for merging:
1. Both ends of one Flye contig aligns to the ends of two different HiFi contigs.
2. The alignment overlap should be minimum 10kb.
3. The orientation is consistent between all three contigs.
4. There is at least one long read spanning the inferred junctions.

Create a new contig by concatenating the two HiFi contigs, using the Flye contig to fill in the gap between them. The coordinates from the alignment are used to extract the correct Flye sequence.
```
samtools faidx sample.hifi.fasta hifi_contig1 > part1.fasta
samtools faidx sample.hifi.fasta hifi_contig2 > part3.fasta
samtools faidx sample.flye.fasta flye_contig:start-end > part2.fasta
cat part1.fasta part2.fasta part3.fasta > merged_contig.fasta
```
`start` and `end` are extracted from the alignment file.
