"""
Creates the dataset containing n-peptide sequences and their features on the pafig dataset
"""
import os
import random

aaimap = {"A":0, "R":1, "N":2, "D":3, "C":4, "Q":5, "E":6, "G":7, 
          "H":8, "I":9, "L":10, "K":11, "M":12, "F":13, "P":14, "S":15, 
          "T":16, "W":17, "Y":18, "V":19} #key-index pairs

def create_pafig_test(n):
    #if os.path.exists("data/test/pafig_hexpepset.txt"):
     #   return
    if not os.path.exists("data/temp/amino_acid_index.txt"):
        create_aaindex()
    f = open("data/test/pafig_dataset.txt")
    data = []
    for line in f:
		if line.strip()[0]=="+":
			data.append(line.split()[1] + " 1") # Positive data
		else:
			data.append(line.split()[1] + " 0") # Negative data
    # Shuffle the data randomly so that we can do cross-validation
    random.shuffle(data)
    amylnset = open("data/test/pafig_hexpepset.txt", "w")
    # Copy the amino acid index to memory and remove incomplete entries
    aai = open("data/temp/amino_acid_index.txt")
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
                fsum += feature[aaimap[x]]
            seq_features = seq_features + " " + str(fsum)
        # sequence info, scaled index value in order
        ac = " "
        for x in data[i].split()[0]:
            k = float(aaimap[x])/19.0*4.0-2.0
            ac = ac + str(k) + " "
        amylnset.write(data[i] + seq_features + ac + "\n")
    print "The pafig_hexpepset.txt has been created."
