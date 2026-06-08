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

