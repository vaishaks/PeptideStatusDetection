"""
Creates the dataset containing n-peptide sequences and their features.
"""

import os
import string
import random
import pandas as pd
import numpy as np

aaimap = {"A":0, "R":1, "N":2, "D":3, "C":4, "Q":5, "E":6, "G":7, 
          "H":8, "I":9, "L":10, "K":11, "M":12, "F":13, "P":14, "S":15, 
          "T":16, "W":17, "Y":18, "V":19} #key-index pairs

def create_aaindex():
    aaindex = open("data/aaindex1.txt")
    amino_acid_index = open("data/temp/amino_acid_index.txt", "w")
    line  = aaindex.readline()
    d = False
    feature = []
    while line:        
        if line[0] == 'D':
	   feature.append("-".join(line[:string.find(line, '(')].split()[1:]))
	   d = True	
        if line[0] == 'I':
            if d == False:
                feature.append("noname")
            f = aaindex.readline().rstrip() + " " + aaindex.readline().rstrip()
            if 'NA' not in f:
                p = map(float, f.split())
                pdash = []
                for v in p:
                    # min-max Normalization, new limts (-2, 2)
                    pdash.append((v - min(p))/(max(p)-min(p))*4 - 2)
                feature.append(" ".join(map(str, pdash)) + "\n")
                amino_acid_index.write(" ".join(feature))  
                feature = []  
            else:
                feature = []            
            d = False
        line = aaindex.readline()
    
    # Adding fifty extra features
    fiftyprops = open("data/50extraprops.txt", "r")
    fiftypropsnames = open("data/50extrapropsnames.txt", "r")
    fiftyprops.readline()
    fiftyprops.readline()
    fiftyprops.readline()

    extrapropmap = [0, 14, 11, 1, 2, 13, 3, 5, 6, 7, 9, 8, 10, 4, 12, 15, 16, 18, 19, 17]

    for (props,name) in zip(fiftyprops, fiftypropsnames):
        props = map(float, props.strip('\r\n').split()[2:])
        name = '-'.join(name.strip('\r\n').split()[2:])
        writeprops = range(20)
        c = 0
        for i in extrapropmap:
            writeprops[c] = (props[i] - min(props))/(max(props) - min(props))*4 - 2
            c += 1  
        amino_acid_index.write(name + " " + " ".join(map(str, writeprops)) + '\n')  

    fiftyprops.close()
    fiftypropsnames.close()
    amino_acid_index.close()
    aaindex.close()
    
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
    f.close()
    
def getPositiveAmylpredData(id, l, r, n):
    """
    Saves all the positive regions for amylpred in a new file 
    if their length is greater than or equal to n.
    """
    if r-l+1 < n:
        return
    am_p = open("data/temp/amylpred_positive_data.txt", "a")
    am_prot = open("data/Positive_amylpred/"+id+".fasta")
    am_prot.readline()
    seq = ''
    for line in am_prot:
        seq += line.rstrip()
    am_p.write(seq[l-1:r]+"\n")
    am_p.close()
    am_prot.close()
    

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
        i += 1

def getData(n):
    if not os.path.exists("data/temp/positive_data.txt"):
        pos = open("data/positive_regions.txt")
        for line in pos:
            l = line.split()
            getPositiveData(l[0], int(l[1]), int(l[2]), n)
    if not os.path.exists("data/temp/amylpred_positive_data.txt"):
        pos = open("data/amylpred_positive_regions.txt")
        for line in pos:
            l = line.split()
            getPositiveAmylpredData(l[0], int(l[1]), int(l[2]), n)            
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
            i += 1
            
def create_amylpred_npeptide_data(n):
    """
    Extracts n-mers from the positive amylpred data and saves them in a new file.
    """
    getData(n)
    f = open("data/temp/amylpred_positive_data.txt")
    npep = open("data/temp/amylpred"+str(n)+"peptides.txt", "w")
    for line in f:
        i = 0
        while i+n < len(line):    
            npep.write(line[i:i+n]+"\n")
            i += 1        
            
def compute_features(seq, feature_ids=[]):
    aai = open("data/temp/amino_acid_index.txt")
    seq_features = []
    aaindex = []
    for line in aai:
        try:
            aaindex.append(map(float, line.split()[1:]))
        except ValueError:
            pass
    for feature in aaindex:
        fsum = 0
        for x in seq:
            fsum += feature[aaimap[x]]
        seq_features.append(fsum)
    for x in seq:
        k = float(aaimap[x])/19.0*4.0-2.0
        seq_features.append(k)
    # Compute all features if feature ID's are not mentioned
    if feature_ids == []:
        return seq_features
    # Using only best 100 features
    seq_features = np.array(seq_features)
    seq_features_trans = seq_features.T
    c = 0
    seq_features = list()
    for i in feature_ids:
        seq_features.append(seq_features_trans[i])
        c += 1
        if c == 100:
            break
    seq_features = np.array(seq_features).T
    return seq_features.tolist()
  
def create_amylnset(n):
    if os.path.exists("data/temp/amyl"+str(n)+"set.txt"):
        return
    if not os.path.exists("data/temp/amino_acid_index.txt"):
        create_aaindex()
    
    create_npeptide_data(n)
    fp = open("data/temp/"+str(n)+"peptides.txt")
    fn = open("data/temp/neg-"+str(n)+"peptides.txt")
    data = [line.rstrip() + " 1" for line in fp.readlines()] # Positive data
    neg = [line.rstrip() + " 0" for line in fn.readlines()] # Negative data
    data.extend(neg)
    # Shuffle the data randomly so that we can do cross-validation
    random.shuffle(data)

    # Creating dataset with all features.
    temp_amylnset = open("data/temp/temp_amyl"+str(n)+"set.txt", "w")
    # Compute the features for each sequence and append them to the data
    for i in xrange(len(data)):        
        seq_features = " ".join(str(e) for e in compute_features(data[i].split()[0]))
        temp_amylnset.write(data[i] + " " + seq_features + "\n")
    temp_amylnset.close()
    
    # Create the .csv file of the sorted scores of features if it does not exist
    if not os.path.exists("data/temp/amylnset_feature_dataframe.csv"):
        import optimal_feature_selection as ofs
        ofs.select_optimal_features(6, "amylnset")
        
    # Creating dataset with optimal features.
    amylnset = open("data/temp/amyl"+str(n)+"set.txt", "w")
    feature_dataframe = pd.read_csv("data/temp/amylnset_feature_dataframe.csv", 
                                        index_col=0, header=0)
    feature_ids = [x for x in feature_dataframe["id"]]
    feature_ids.extend(range(len(feature_ids), len(feature_ids)+n))
        
    # Compute the features for each sequence and append them to the data
    for i in xrange(len(data)):        
        seq_features = " ".join(str(e) for e in compute_features(data[i].split()[0],
                                                                    feature_ids))
        amylnset.write(data[i] + " " + seq_features + "\n")

def create_pafig_data(n):
    if os.path.exists("data/temp/pafig_hexpepset.txt"):
        print "Using existing data."
        return
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
    
    # Creating a dataset with all features
    temp_pafig_hexpepset = open("data/temp/temp_pafig_hexpepset.txt", "w")   
    # Compute the features for each sequence and append them to the data
    for i in xrange(len(data)):        
        seq_features = " ".join(str(e) for e in compute_features(data[i].split()[0]))
        temp_pafig_hexpepset.write(data[i] + " " + seq_features + "\n")
    temp_pafig_hexpepset.close()
    
    if not os.path.exists("data/temp/pafig_feature_dataframe.csv"):
        import optimal_feature_selection as ofs
        ofs.select_optimal_features(6, "pafig")
    
    feature_dataframe = pd.read_csv("data/temp/pafig_feature_dataframe.csv", 
                                        index_col=0, header=0)
    feature_ids = [x for x in feature_dataframe["id"]]
    feature_ids.extend(range(len(feature_ids), len(feature_ids)+n))       
    
    # Compute the features for each sequence and append them to the data
    pafig_hexpepset = open("data/temp/pafig_hexpepset.txt", "w") 
    for i in xrange(len(data)):        
        seq_features = " ".join(str(e) for e in compute_features(data[i].split()[0],
                                                                    feature_ids))
        pafig_hexpepset.write(data[i] + " " + seq_features + "\n")
    print "The pafig_hexpepset.txt has been created."

def create_zipper_data(n):
    if os.path.exists("data/temp/zipper_hexpepset.txt"):
        print "Using existing data."
        return
    if not os.path.exists("data/temp/amino_acid_index.txt"):
        create_aaindex()
    
    f = open("data/test/zipper_dataset.txt")
    data = []
    for line in f:
        if line.strip()[0]=="+":
            data.append(line.split()[1] + " 1") # Positive data
        else:
            data.append(line.split()[1] + " 0") # Negative data
    # Shuffle the data randomly so that we can do cross-validation
    random.shuffle(data)
    
    # Creating a dataset with all features
    temp_zipper_hexpepset = open("data/temp/temp_zipper_hexpepset.txt", "w")   
    # Compute the features for each sequence and append them to the data
    for i in xrange(len(data)):        
        seq_features = " ".join(str(e) for e in compute_features(data[i].split()[0]))
        temp_zipper_hexpepset.write(data[i] + " " + seq_features + "\n")
    temp_zipper_hexpepset.close()
    
    if not os.path.exists("data/temp/zipper_feature_dataframe.csv"):
        import optimal_feature_selection as ofs
        ofs.select_optimal_features(6, "zipper")
    
    feature_dataframe = pd.read_csv("data/temp/zipper_feature_dataframe.csv", 
                                        index_col=0, header=0)
    feature_ids = [x for x in feature_dataframe["id"]]
    feature_ids.extend(range(len(feature_ids), len(feature_ids)+n))  
    
    # Compute the features for each sequence and append them to the data
    zipper_hexpepset = open("data/temp/zipper_hexpepset.txt", "w") 
    for i in xrange(len(data)):        
        seq_features = " ".join(str(e) for e in compute_features(data[i].split()[0],
                                                                    feature_ids))
        zipper_hexpepset.write(data[i] + " " + seq_features + "\n")
    print "The zipper_hexpepset.txt has been created."

def create_amylpred_data(n):
    if os.path.exists("data/temp/amylpred"+str(n)+"set.txt"):
        print "Using existing data."
        return
    if not os.path.exists("data/temp/amino_acid_index.txt"):
        create_aaindex()
    
    create_amylpred_npeptide_data(n)
    fp = open("data/temp/amylpred"+str(n)+"peptides.txt")
    fn = open("data/temp/neg-"+str(n)+"peptides.txt")
    data = [line.rstrip() + " 1" for line in fp.readlines()] # Positive data
    neg = [line.rstrip() + " 0" for line in fn.readlines()] # Negative data
    data.extend(neg)
    # Shuffle the data randomly so that we can do cross-validation
    random.shuffle(data)
    
    # Creating a dataset with all features
    temp_amylprednset = open("data/temp/temp_amylpred"+str(n)+"set.txt", "w")   
    # Compute the features for each sequence and append them to the data
    for i in xrange(len(data)):        
        seq_features = " ".join(str(e) for e in compute_features(data[i].split()[0]))
        temp_amylprednset.write(data[i] + " " + seq_features + "\n")
    temp_amylprednset.close()
    
    if not os.path.exists("data/temp/amylpred_feature_dataframe.csv"):
        import optimal_feature_selection as ofs
        ofs.select_optimal_features(n, "amylpred")
    
    feature_dataframe = pd.read_csv("data/temp/amylpred_feature_dataframe.csv", 
                                        index_col=0, header=0)
    feature_ids = [x for x in feature_dataframe["id"]]
    feature_ids.extend(range(len(feature_ids), len(feature_ids)+n))  
    
    # Compute the features for each sequence and append them to the data
    amylprednset = open("data/temp/amylpred"+str(n)+"set.txt", "w") 
    for i in xrange(len(data)):        
        seq_features = " ".join(str(e) for e in compute_features(data[i].split()[0],
                                                                    feature_ids))
        amylprednset.write(data[i] + " " + seq_features + "\n")
    print "The amylprednset.txt has been created."
