import sys
from utils import *

# extract genome seq from vcf 

vcffile = sys.argv[1]
prot_table = sys.argv[2]

sites = {}
for line in open(vcffile):
  if line[0]=='#': continue
  w = line.rstrip().split('\t')
  pos,ref,alt,code = int(w[1]),w[3],w[4],w[6] # 1-based
  if alt=="<DUP>": continue # just skip it, coord defined on previous line
  sites[pos] = (ref[0],alt[0],code) # if multiple nucs (ins or del), just keep first nucs for ref and alt

for line in open(prot_table):
  w = line.rstrip().split('\t')
  rv,gene,start,end,strand = w[8],w[7],int(w[1]),int(w[2])+1,w[3] # coords are 1-based; make end exclusive

  s = ""
  for i in range(start,end):
    ref,alt,code = sites[i]
    if "Del" in code: nuc = '-'
    #elif "LowCov" in code: nuc = '?'
    elif "LowCov" in code: nuc = ref.lower() # instead of '?', ignore SNPs and Amb
    elif "PASS" not in code: nuc = '#' # mostly Amb
    elif alt=='.': nuc = ref
    else: nuc = alt # this works for indels too, since just kept first nuc above
    s += nuc
  if strand=='-': s = reverse_complement(s)

  print("> %s" % rv)
  print(s)
 
  expected_len = end-start
  actual_len = len(s)
  delta = actual_len-expected_len
  vals = [vcffile,rv,gene,expected_len,actual_len,delta]
  #sys.stderr.write("%s\n" % ('\t'.join([str(x) for x in vals])))
