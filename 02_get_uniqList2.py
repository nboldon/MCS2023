import sys
import gzip
file1 = "/home/qisun/sc_atac/737K-cratac-v1.txt.gz"
file2 = "/home/qisun/sc_atac/shuf_737K-cratac-v1.txt.gz"
def hamming_distance(string1, string2):
    # Start with a distance of zero, and count up
    distance = 0
    # Loop over the indices of the string
    L = len(string1)
    for i in range(L):
        # Add 1 to the distance if these two characters are not equal
        if string1[i] != string2[i]:
            distance += 1
    # Return the final count of differences
    return distance


whiteBarcodes = []
with gzip.open(file1, "rb") as IN:
    for line in IN:
        whiteBarcodes.append(line.decode('utf-8').strip())
IN.close()

skip = 40020
max = 500000
i=0
GGG =open("extra", "wt")
with gzip.open(file2, "rb") as IN:
    for line in IN:
        i+=1
        if i>skip:
            break
    for line in IN:
        i+=1
        newbc = line.decode('utf-8').strip()
        if i>max:
            break
        flag=True    
        for t in whiteBarcodes:
            if(hamming_distance(t, newbc) <3):
                flag = False
                break
        if (flag):
            GGG.write(newbc+"\n")
        
IN.close()
GGG.close()


