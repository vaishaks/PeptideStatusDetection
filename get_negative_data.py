"""
Creates the negative n-peptide data.
"""

import os

def getNegativeData(id, n):
    npep = open("data/temp/neg-"+str(n)+"peptides.txt", "a")
    prot = open("data/Negative/"+id+".fasta")
    prot.readline()
    for line in prot:
        seq = ''
        for line in prot:
            seq += line.rstrip()
    i = 0
    while i+n < len(seq):    
        npep.write(seq[i:i+6]+"\n")
        i += n

def start(n):
    if not os.path.exists("data/temp/neg-"+str(n)+"peptides.txt"):
        f = open("data/negative_regions.txt")
        for line in f:
            getNegativeData(line.rstrip(), n)
            