"""
Extracts n-mers from the positive data and saves them in a new file.
"""

import get_positive_data as gp
import get_negative_data as gn

def create_npeptide_data(n):
    gp.start(n)
    gn.start(n)
    f = open("data/temp/positive_data.txt")
    npep = open("data/temp/"+str(n)+"peptides.txt", "w")
    for line in f:
        i = 0
        while i+n < len(line):    
            npep.write(line[i:i+6]+"\n")
            i += n

create_npeptide_data(7)