import sys
import gzip
import os

sample = 'Landrace'

chromlen = {}
infile = f'{sample}.chromlen.tsv'
# Read the length of each chromosome
with open(infile, 'r') as fin:
    for line in fin:
        l = line.strip().split()
        length = int(l[1])
        chromlen[l[0]] = [1,length]

satellite = {'AJ491340.1':'repeat',
             'AJ889249.1':'Ssc14','AJ920048.1':'Ssc14','AJ920049.1':'Ssc14','AJ920050.1':'Ssc14','AJ920051.1':'Ssc14','AJ920052.1':'Ssc14',
             'AJ937272.1':'Ssc14','AJ937273.1':'Ssc14','AJ937274.1':'Ssc14','AJ937275.1':'Ssc14','AJ937276.1':'Ssc14',
             'U42362.1':'Mc2','U42363.1':'Mc2','U42364.1':'Mc2','U42365.1':'Mc2','U42366.1':'Mc2','U42367.1':'Mc2','U42368.1':'Mc2','U42369.1':'Mc2','U42370.1':'Mc2',
             'X16513.1':'acrocentric',
             'X51555.1':'Mc1','X51556.1':'Mc1','X51557.1':'Mc1','X51558.1':'Mc1','X51559.1':'Mc1','X51560.1':'Mc1',
             'X51561.1':'Ac2','X51562.1':'Ac2','X51563.1':'Ac2','X51564.1':'Ac2','X51565.1':'Ac2',
             'X70941.1':'acrocentric','X62140.1':'metacentric','X62139.1':'metacentric','X62138.1':'metacentric',
             'D11085.1':'acrocentric'}
satellites = ['AJ491340.1','AJ889249.1','AJ920048.1','AJ920049.1','AJ920050.1','AJ920051.1','AJ920052.1',
             'AJ937272.1','AJ937273.1','AJ937274.1','AJ937275.1','AJ937276.1',
             'U42362.1','U42363.1','U42364.1','U42365.1','U42366.1','U42367.1','U42368.1','U42369.1','U42370.1',
             'X16513.1',
             'X51555.1','X51556.1','X51557.1','X51558.1','X51559.1','X51560.1',
             'X51561.1','X51562.1','X51563.1','X51564.1','X51565.1',
             'X70941.1','X62140.1','X62139.1','X62138.1',
             'D11085.1']

# Read the alignment file and record each base covered by at least one alignment
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

# Count number of aligned bases per 10 kbp
step = 10000
CHROM = [] # Chromosome
POS = [] # Mid point of window
REP = [] # Fraction of mapped bases within the window
SAT = [] # Satellite
TYP = [] # Type of satellite
for index,sat in enumerate(satellites):
    infile = f'{sample}/lastz_{sample}_satellites.{sat}.map.gz'
    with gzip.open(infile, 'rt') as fin:
        for line in fin:
            chrom, seq = line.strip().split()
            x1 = 0
            x2 = x1+step
            while x1 < len(seq):
                segment = seq[x1:x2]
                value = segment.count('1') / step
                SAT.append(sat)
                CHROM.append(chrom)
                POS.append((x1+x2)//2)
                REP.append(value)
                TYP.append(satellite[sat])
                x1 = x2
                x2 = x1+step
df = pd.DataFrame({'SAT':SAT,'CHROM':CHROM,'POS':POS,'REP':REP,'TYPE':TYP})
df.to_csv(f'lastz_{sample}_satellites.dataframe-{step}.tsv', sep='\t', index=False)
