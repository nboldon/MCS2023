import sys
import gzip
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

whiteList = "/home/qisun/sc_atac/737K-cratac-v1.txt.gz"
whiteBarcodes = []
max=40000
i=0
with gzip.open(whiteList, "rb") as IN:
    for line in IN:
        i+=1
        whiteBarcodes.append(line.decode('utf-8').strip())
        if i>max:
            break
IN.close()

GGG = open ("goodlist", "wt")
for i in range(0, max-1):
    ccc =0 
    for j in range(i+1, max):
        h =hamming_distance(whiteBarcodes[i], whiteBarcodes[j])
        if h<3:
            ccc +=1
    if ccc==0:
        GGG.write(whiteBarcodes[i]+"\n")
GGG.close()