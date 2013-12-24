"""
Creates the dataset containing n-peptide sequences and their features.
"""

import os
import random

aaimap = {"A":0, "R":1, "N":2, "D":3, "C":4, "Q":5, "E":6, "G":7, 
          "H":8, "I":9, "L":10, "K":11, "M":12, "F":13, "P":14, "S":15, 
          "T":16, "W":17, "Y":18, "V":19} #key-index pairs

def create_aaindex():
    aaindex = open("data/aaindex1.txt")
    amino_acid_index = open("data/amino_acid_index.txt", "w")
    line  = aaindex.readline()
    while line:
        if line[0] == 'I':
            amino_acid_index.write(aaindex.readline().rstrip()
            + " " + aaindex.readline())
            aaindex.readline()
        line = aaindex.readline()
        
def getPositiveData(id, l, r, n):
    """
    Saves all the positive regions in the proteins in a new file 
    if their length is greater than or equal to n.
    """
    if r-l+1 < n:
        return
    f = open("data/temp/positive_data.txt", "a")
    prot = open("data/Positive/"+id+".fasta")
    prot.readline()
    seq = ''
    for line in prot:
        seq += line.rstrip()
    f.write(seq[l-1:r]+"\n")

def getNegativeData(id, n):
    """
    Creates the negative n-peptide data.
    """
    npep = open("data/temp/neg-"+str(n)+"peptides.txt", "a")
    prot = open("data/Negative/"+id+".fasta")
    prot.readline()
    for line in prot:
        seq = ''
        for line in prot:
            seq += line.rstrip()
    i = 0
    while i+n < len(seq):    
        npep.write(seq[i:i+n]+"\n")
        i += n

def getData(n):
    if not os.path.exists("data/temp/positive_data.txt"):
        pos = open("data/positive_regions.txt")
        for line in pos:
            l = line.split()
            getPositiveData(l[0], int(l[1]), int(l[2]), n)
    if not os.path.exists("data/temp/neg-"+str(n)+"peptides.txt"):
        f = open("data/negative_regions.txt")
        for line in f:
            getNegativeData(line.rstrip(), n)
            
def create_npeptide_data(n):
    """
    Extracts n-mers from the positive data and saves them in a new file.
    """
    getData(n)
    f = open("data/temp/positive_data.txt")
    npep = open("data/temp/"+str(n)+"peptides.txt", "w")
    for line in f:
        i = 0
        while i+n < len(line):    
            npep.write(line[i:i+n]+"\n")
            i += n
            
def compute_features(seq):
    aai = open("data/amino_acid_index.txt")
    seq_features = []
    aaindex = []
    for line in aai:
        if len(line.split()) == 20:
            try:
                aaindex.append(map(float, line.split()))
            except ValueError:
                pass
    for feature in aaindex:
        fsum = 0
        for x in seq:
            fsum += feature[aaimap[x]]
        seq_features.append(fsum)
    return seq_features

def create_amylnset(n):
    if os.path.exists("data/temp/amyl"+str(n)+"set.txt"):
        return
    if not os.path.exists("data/amino_acid_index.txt"):
        create_aaindex()
    create_npeptide_data(n)
    fp = open("data/temp/"+str(n)+"peptides.txt")
    fn = open("data/temp/neg-"+str(n)+"peptides.txt")
    data = [line.rstrip() + " 1" for line in fp.readlines()] # Positive data
    neg = [line.rstrip() + " 0" for line in fn.readlines()] # Negative data
    data.extend(neg)
    # Shuffle the data randomly so that we can do cross-validation
    random.shuffle(data)
    amylnset = open("data/temp/amyl"+str(n)+"set.txt", "w")
    # Copy the amino acid index to memory and remove incomplete entries
    aai = open("data/amino_acid_index.txt")
    aaindex = []
    for line in aai:
        if len(line.split()) == 20:
            try:
                aaindex.append(map(float, line.split()))
            except ValueError:
                pass
    # Compute the features for each sequence and append them to the data
    for i in xrange(len(data)):
        seq_features = ""
        for feature in aaindex:
            fsum = 0
            for x in data[i].split()[0]:
                fsum += feature[aaimap[x]] # summation, no auto-correlation
            seq_features = seq_features + " " + str(fsum)
        amylnset.write(data[i] + seq_features + "\n")

create_amylnset(7)