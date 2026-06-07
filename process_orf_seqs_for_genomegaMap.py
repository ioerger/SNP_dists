import sys
sys.path.append("/home/ioerger/genomics")
from utils import *

seqsfile = sys.argv[1]
tag = sys.argv[2]

H,S = read_fasta(seqsfile)

g = read_genome("/mnt/raid/Moldova/H37Rv3.fna")
for line in open("/mnt/raid/Moldova/H37Rv3.prot_table"):
  if tag in line:
    w = line[:-1].split("\t")
    start,end,strand = int(w[1])-1,int(w[2])-1,w[3]
    orf,gene = w[8],w[7]
    refseq = g[start:end+1]
    if strand=='-': refseq = reverse_complement(refseq)
    refseq = refseq[:-3] # trim-off stop codon

STOP_CODONS = "TAG TGA TAA".split()

for h,s in zip(H,S):
  s = s[:-3] # note, assume input has a stop codon
  #s = s.upper() # ignore LowCov
  # note: by leaving as lower-case, will get converted to '-' below
  #if s[-3:] in STOP_CODONS: s = s[:-3].upper()
  i,n = 0,len(s)
  s2 = ""
  while i<n:
    triplet = s[i:i+3]
    aa = codon.get(triplet,'X')
    if aa in 'X*': triplet = '---' # ambig or stop codon
    s2 += triplet
    i += 3

  # remove SNP clusters (>=3 SNPs within 50bp, relative to H37Rv3)
  if False:
    W,C = 50,2 # W is window size in bp, C is max SNPs within window allowed
    SNPs = []
    for i in range(len(s2)):
      if s2[i]!=refseq[i] and s2!='-' and s2!='#': SNPs.append(i)
    nSNPs = len(SNPs)
    badSNPs = {}
    for p in range(nSNPs-C):
      q = p+C
      if SNPs[q]-SNPs[p]<W:
        for r in range(p,q+1): badSNPs[SNPs[r]] = 1
    badSNPs = sorted(badSNPs.keys())
    if len(badSNPs)>0: sys.stderr.write("removing SNP cluster in %s in %s (%s in %s-%s)\n" % (tag,h[1:],len(badSNPs),badSNPs[0],badSNPs[-1]))
    for r in badSNPs: s2 = s2[:r]+refseq[r]+s2[r+1:]

  print(h)
  print(s2)
