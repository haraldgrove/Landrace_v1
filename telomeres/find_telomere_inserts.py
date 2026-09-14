import sys
import gzip

sample = 'Landrace'

outfile = f'{sample}.inbetweens.front.tsv'

infile = f'{sample}.chromlen.tsv'
chroms = []
with open(infile, 'r') as fin:
    for line in fin:
        chrom,length = line.strip().split()
        chroms.append(chrom)

pattern = 'CCCTAA'

with open(outfile, 'w') as fout:
    for index, chrom in enumerate(chroms):
        infile = f'{chrom}.telomere.fasta'
        with open(infile, 'r') as fin:
            header = next(fin)
            seq1 = next(fin).strip().upper()
            pos = 0
            extra = ''
            while pos < len(seq1):
                pat = seq1[pos:pos+6] # current 6-base pattern to consider
                if pat == pattern:
                    if len(extra) > 0:
                        fout.write(f'{index+1}\t{pos-len(extra)}\t{pos}\t{pattern}\t{extra}\n')
                        extra = ''
                    pos += 6
                    continue
                extra += seq1[pos]
                pos += 1
            else:
                if len(extra) > 0:
                    fout.write(f'{index+1}\t{pos-len(extra)}\t{pos}\t{pattern}\t{extra}\n')

def reverse_complement(dna_seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N':'N'}
    reversed_seq = dna_seq[::-1]
    reverse_complement_seq = ''.join([complement[base] for base in reversed_seq])
    return reverse_complement_seq

outfile = f'{sample}.inbetweens.rev_comp.end.tsv'

pattern = 'CCCTAA'

with open(outfile, 'w') as fout:
    for index, chrom in enumerate(chroms):
        infile = f'{chrom}.telomere.fasta'
        with open(infile, 'r') as fin:
            header = next(fin)
            seq1 = next(fin)
            seq2 = reverse_complement(next(fin).strip().upper())
            pos = 0
            extra = ''
            while pos < len(seq2):
                pat = seq2[pos:pos+6] # current 6-base pattern to consider
                if pat == pattern:
                    if len(extra) > 0:
                        fout.write(f'{index+1}\t{pos-len(extra)}\t{pos}\t{pattern}\t{extra}\n')
                        extra = ''
                    pos += 6
                    continue
                extra += seq2[pos]
                pos += 1
            else:
                if len(extra) > 0:
                    fout.write(f'{index+1}\t{pos-len(extra)}\t{pos}\t{pattern}\t{extra}\n')
