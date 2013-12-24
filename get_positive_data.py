"""
Saves all the positive regions in the proteins in a new file if their length is 
greater than or equal to n.
"""

import os

def getPositiveData(id, l, r, n):
    if r-l+1 < n:
        return
    f = open("data/temp/positive_data.txt", "a")
    prot = open("data/Positive/"+id+".fasta")
    prot.readline()
    seq = ''
    for line in prot:
        seq += line.rstrip()
    f.write(seq[l-1:r]+"\n")

def start(n):
    if not os.path.exists("data/temp/positive_data.txt"):
        pos = open("data/positive_regions.txt")
        for line in pos:
            l = line.split()
            getPositiveData(l[0], int(l[1]), int(l[2]), n)
