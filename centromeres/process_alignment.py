import sys
import gzip
import os

sample = 'Landrace'

try:
    os.mkdir(sample)
except FileExistsError:
    pass
     
chromlen = {}
infile = f'{sample}.chromlen.tsv'
with open(infile, 'r') as fin:
    for line in fin:
        l = line.strip().split()
        length = int(l[1])
        chromlen[l[0]] = [1,length]

satellite = {'AJ491340.1':'repeat',
             'AJ889249.1':'Ssc14','AJ920048.1':'Ssc14','AJ920049.1':'Ssc14','AJ920050.1':'Ssc14','AJ920051.1':'Ssc14','AJ920052.1':'Ssc14','AJ937272.1':'Ssc14','AJ937273.1':'Ssc14','AJ937274.1':'Ssc14','AJ937275.1':'Ssc14','AJ937276>
             'U42363.1':'Mc2','U42364.1':'Mc2','U42365.1':'Mc2','U42366.1':'Mc2','U42367.1':'Mc2','U42368.1':'Mc2','U42369.1':'Mc2','U42370.1':'Mc2',
             'X16513.1':'repeat',
             'X51555.1':'Mc1','X51556.1':'Mc1','X51557.1':'Mc1','X51558.1':'Mc1','X51559.1':'Mc1','X51560.1':'Mc1',
             'X51561.1':'Ac2','X51562.1':'Ac2','X51563.1':'Ac2','X51564.1':'Ac2','X51565.1':'Ac2',
             'X70941.1':'centromere',
             'D11085.1':'acrocentric',
             'X62140.1':'centromere',
             'X62139.1':'centromere',
             'X62138.1':'centromere'}

infile = f'{sample}_satellites.paf.gz'
for query in satellite:
    outfile = f'{sample}_satellites.{query}.map.gz'
    repeats = {}
    for chrom in chromlen:
        x1,x2 = chromlen[chrom]
        repeats[chrom] = ['0']*(x2-x1+1)
    with gzip.open(infile, 'rt') as fin:
        header = next(fin)
        for line in fin:
            qname,qlen,qstart,qend,strand,tname,tlen,tstart,tend,match,length,mapq,*tags = line.strip().split()
            if tname not in chromlen:
                continue
            if qname != query:
                continue
            seqid = int(match)/int(length)
            coverage = (int(qend)-int(qstart))/int(qlen)
            if coverage < 0.9 or seqid < 0.6:
                continue
            START = int(tstart)
            END = int(tend)
            for i in range(START,END+1):
                repeats[tname][i-1] = '1'
    with gzip.open(outfile, 'wt') as fout:
        for chrom in chromlen:
            seqlist = repeats[chrom]
            outseq = ''.join(seqlist)
            fout.write(f'{chrom}\t{outseq}\n')
