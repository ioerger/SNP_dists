import sys

seqsfile = sys.argv[1]
tag = sys.argv[2]

###################

def read_fasta(filename):
  headers,seqs = [],[]
  seq = ""
  for line in open(filename):
    line = line.rstrip()
    if len(line)==0: continue
    if line[0]==">":
      headers.append(line)
      if seq!="": seqs.append(seq)
      seq = ""
    else: seq += line
  seqs.append(seq)
  return headers,seqs

def read_genome(fname):
  h,s = read_fasta(fname)
  return s[0]

complement = {'A':'T','T':'A','C':'G','G':'C'}

def reverse_complement(seq):
  s = list(seq)
  s.reverse()
  for i in range(len(s)):
    s[i] = complement.get(s[i],s[i]) # if unknown, leave as it, e.g > or !
  s = ''.join(s)
  return s

def translate(seq):
  s = ""
  for i in range(int(len(seq)/3)):
    s += codon.get(seq[3*i:3*(i+1)],'X')
  return s

codon    = { 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
             'CTA': 'L', 'CTG': 'L', 'CTN': 'L', 'TGG': 'W',
             'TAA': '*', 'TAG': '*', 'TGA': '*', 'ATG': 'M',
             'TTT': 'F', 'TTC': 'F', 'TAT': 'Y', 'TAC': 'Y',
             'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
             'TCN': 'S', 'AGT': 'S', 'AGC': 'S', 'CCT': 'P',
             'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CCN': 'P',
             'TGT': 'C', 'TGC': 'C', 'CAT': 'H', 'CAC': 'H',
             'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
             'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
             'CGN': 'R', 'AGA': 'R', 'AGG': 'R', 'ATT': 'I',
             'ATC': 'I', 'ATA': 'I', 'AAA': 'K', 'AAG': 'K',
             'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
             'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
             'ACN': 'T', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
             'GTG': 'V', 'GTN': 'V', 'GCT': 'A', 'GCC': 'A',
             'GCA': 'A', 'GCG': 'A', 'GCN': 'A', 'GGT': 'G',
             'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GGN': 'G',
             'TAN': 'X', 'TTN': 'X', 'TGN': 'X', 'CAN': 'X',
             'ATN': 'X', 'AAN': 'X', 'GAN': 'X', 'AGN': 'X',
             'ANA': 'X', 'ANT': 'X', 'ANG': 'X', 'ANC': 'X',
             'TNA': 'X', 'TNT': 'X', 'TNG': 'X', 'TNC': 'X',
             'GNA': 'X', 'GNT': 'X', 'GNG': 'X', 'GNC': 'X',
             'CNA': 'X', 'CNT': 'X', 'CNG': 'X', 'CNC': 'X',
             'NAA': 'X', 'NAT': 'X', 'NAG': 'X', 'NAC': 'X',
             'NTA': 'X', 'NTT': 'X', 'NTG': 'X', 'NTC': 'X',
             'NGA': 'X', 'NGT': 'X', 'NGG': 'X', 'NGC': 'X',
             'NCA': 'X', 'NCT': 'X', 'NCG': 'X', 'NCC': 'X',
             'NNA': 'X', 'NNT': 'X', 'NNG': 'X', 'NNC': 'X',
             'ANN': 'X', 'TNN': 'X', 'GNN': 'X', 'NNC': 'X',
             'NAN': 'X', 'NTN': 'X', 'NGN': 'X', 'NCN': 'X',
             'NNN': 'X'}

###################

H,S = read_fasta(seqsfile)

g = read_genome("H37Rv3.fna")
for line in open("H37Rv3.prot_table"):
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
