def getPositiveData(id, l, r, n):
    if r-l+1 < n:
        return
    f = open("data/positive_data.txt", "a")
    prot = open("data/Positive/"+id+".fasta")
    prot.readline()
    seq = ''
    for line in prot:
        seq += line.rstrip()
    f.write(seq[l-1:r]+"\n")

def start(n):
    pos = open("data/positive_regions.txt")
    for line in pos:
        l = line.split()
        getPositiveData(l[0], int(l[1]), int(l[2]), n)
