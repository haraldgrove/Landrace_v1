#!/bin/bash

target=Sus_scrofa_gca963921485v1.norwegian_landrace.dna.toplevel.fa.gz
query=sscro_satellites.fasta
chrom=$1
output=${chrom}_satellites.paf

singularity exec lastz1.0.4 \
    lastz \
    ${target}[subset=${chrom}] \
    ${query} \
    --output=${output} \
    --ambiguous=iupac \
    --format=PAF
