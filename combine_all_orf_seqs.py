import sys
from utils import *

# example:
# > python3 ~/SNP_dists/combine_all_orf_seqs.py all.genes strains.batchE
# reads <strain>.orfs.txt; writes <orf>.seqs.txt

all_orfs = sys.argv[1]
all_strains = sys.argv[2]

OrfSeqs = {}

for line in open(all_strains):
  id = line[:-1]
  fname = "%s.orfs.fna" % id
  sys.stderr.write("reading %s\n" % fname)
  h,s = read_fasta(fname)
  for h2,s2 in zip(h,s):
    #print("%s  %s" % (id,s2)) # phylip
    #print(">%s\n%s" % (id,s2)) # fasta
    orf = h2.split()[1]
    if orf not in OrfSeqs: OrfSeqs[orf] = []
    OrfSeqs[orf].append((id,s2))

for line in open(all_orfs):
  orf = line[:-1]
  fname = "%s.seqs.txt" % orf
  sys.stderr.write("writing %s\n" % fname)
  fil = open(fname,"w+")
  for (strain,seq) in OrfSeqs[orf]:
    fil.write(">%s\n%s\n" % (strain,seq)) # fasta  
  fil.close()
