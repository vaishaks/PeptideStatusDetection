aaimap = {"A":0, "R":1, "N":2, "D":3, "C":4, "Q":5, "E":6, "G":7, "H":8, "I":9,
          "L":10, "K":11, "M":12, "F":13, "P":14, "S":15, "T":16, "W":17, "Y":18,
          "V":19} #key-index pairs

def create():
    aaindex = open("data/aaindex1.txt")
    amino_acid_index = open("data/amino_acid_index.txt", "w")
    line  = aaindex.readline()
    while line:
        if line[0] == 'I':
            amino_acid_index.write(aaindex.readline().rstrip()
            + " " + aaindex.readline())
            aaindex.readline()
        line = aaindex.readline()
    