#! /usr/bin/env python

import sys
import gzip

infile = sys.argv[1] # Fasta file with one header line and one sequence line per entry.
outfile = sys.argv[2] # Columns: chrom, start_pos, number of telomere repeat units
pattern1 = 'TTAGGG'
pattern2 = 'CCCTAA'

count = 0
with gzip.open(infile,'rt') as fin, open(outfile, 'w') as fout:
    for line in fin:
        if line.startswith('>'):
            chrom = line.strip()[1:]
            continue
        pos = 0
        seq = line.strip().upper()
        while True:
           start = seq.find(pattern1, pos)
            i = start
            if i == -1:
                break
            count = 1
            while seq[i+6:i+12] == pattern1:
                count += 1
                i += 6
            else:
                fout.write(f'{chrom}\t{start}\t{count}\n')
                pos = i+6
        pos = 0
        seq = line.strip().upper()
        while True:
            start = seq.find(pattern2, pos)
            i = start
            if i == -1:
                break
            count = 1
            while seq[i+6:i+12] == pattern2:
                count += 1
                i += 6
            else:
                fout.write(f'{chrom}\t{start}\t{count}\n')
                pos = i+6
