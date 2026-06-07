import sys

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

###################

all_orfs = sys.argv[1]
all_strains = sys.argv[2]

OrfSeqs = {}

for line in open(all_strains):
  id = line[:-1]
  fname = "%s/%s.orfs.fna" % (id,id)
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
  fname = "orfs/%s.seqs.txt" % orf
  sys.stderr.write("writing %s\n" % fname)
  fil = open(fname,"w+")
  for (strain,seq) in OrfSeqs[orf]:
    fil.write(">%s\n%s\n" % (strain,seq)) # fasta  
  fil.close()
