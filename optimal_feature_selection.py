import os
from sklearn import svm

aaindex = open("data/aaindex1.txt")
amino_acid_index = open("data/temp/amino_acid_index.txt", "w")
line  = aaindex.readline()
while line:
	if line[0] == 'I':
		f = aaindex.readline().rstrip() + " " + aaindex.readline().rstrip()
		if 'NA' not in f:
			p = map(float, f.split())
			pdash = []
			for v in p:
				# min-max Normalization, new limts (-2, 2)
				pdash.append((v - min(p))/(max(p)-min(p))*4 - 2)
			amino_acid_index.write(" ".join(map(str, pdash)))
		aaindex.readline()
	if line[0] == 'D':
		f = " ".join(line.split()[1:])
		amino_acid_index.write(" "+ f + "\n")
	line = aaindex.readline()
