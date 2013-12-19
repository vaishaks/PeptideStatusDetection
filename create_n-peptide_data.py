import getPositiveData as g

def create_npeptide_data(n):
    g.start(n)
    f = open("data/positive_data.txt")
    npep = open("data/"+str(n)+"peptides.txt", "w+")
    for line in f:
        i = 0
        while i+6 < len(line):    
            npep.write(line[i:i+6])

create_npeptide_data(6)