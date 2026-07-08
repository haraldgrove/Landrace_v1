# Landrace_v1
Scripts used for the creation and analysis of the Landrace_v1 assembly.

## Assembly

### Flye

Reads were filtered using Filtlong:
```
filtlong --min_length 4000 --keep_percent 90 sample.fq.gz | pigz -p4 > sample.filtered.fq.gz
```

Filtered ONT reads were assembled into contigs using Flye (v 2.9.1-b1780) 
```
flye \
  --nano-raw sample.filtered.fq.gz \
  --min-overlap 20000 \
  --genome-size 2.5g \
  --out-dir sample.flye.outdir \
  --asm-coverage 40 \
  -t 30
```
The “--asm-coverage” flag was set to 40 to reduce the memory consumption in the initial disjoining step
A range of required read overlaps (--min-overlap flag) were used (5k, 7k, 10k, 15k, 20k, 25k and 30k) to produce multiple independent assemblies. 

### HiFIasm

A different assembly was created with ONT reads and HiFi reads using HiFiasm v 0.19.6-r595 default settings.
```
hifiasm -o sample_hifiasm -t 30 --ul sample.filtered.fq.gz sample.hifi.fq.gz
```

### Merging
Align HiFi assembly contigs against Flye assembly contigs.
```
<software> sample.flye.fasta sample.hifi.fasta
```

Identify possible cases for merging:
1. Both ends of one Flye contig aligns to the ends of two different HiFi contigs.
2. The alignment overlap should be minimum 10kb.
3. The orientation is consistent between all three contigs.
4. There is at least one long read spanning the inferred junctions.

Create a new contig by concatenating the two HiFi contigs, using the Flye contig to fill in the gap between them. The coordinates from the alignment are used to extract the correct Flye sequence.
