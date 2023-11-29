from loading import load_directory
from kmers import stream_kmers, kmer2str
import numpy as np
import time
import sys

code = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

def jaccard(fileA, fileB, k):
    """
    fileA / fileB : lists of lists containing sequences on different chromosomes in A or B
    k : integer for the length of kmers
    """
    kmers_A = {}
    for list_Ai in fileA:
        for kmer,kmer_comp in stream_kmers(list_Ai, k):
            kmers_A[kmer]      = kmers_A.get(kmer, 0) + 1
            kmers_A[kmer_comp] = kmers_A.get(kmer_comp, 0) + 1
            
    kmers_B = {}
    for list_Bi in fileB:
        for kmer,kmer_comp in stream_kmers(list_Bi, k):
            kmers_B[kmer]      = kmers_B.get(kmer, 0) + 1
            kmers_B[kmer_comp] = kmers_B.get(kmer_comp, 0) + 1

    intersect = 0
    for key in kmers_A:
        intersect += min( kmers_A[key], kmers_B.get(key,0) )

    #union = sum(kmers_A.values()) + sum(kmers_B.values()) 
    union = 0
    for key in kmers_A:
        union += max( kmers_A[key], kmers_B.get(key,0) )
    for key in kmers_B:
        union += max( kmers_B[key], kmers_A.get(key,0) )
    # A.U.B = A+B-(A/-\B)
    union -= intersect
    
    size = sys.getsizeof(kmers_A) + sys.getsizeof(kmers_B)
    del kmers_A
    del kmers_B
    
    return intersect / union, size



if __name__ == "__main__":
    # Load all the files in a dictionary
    start = time.time()
    files = load_directory("data")
    k = 21
    l = len(files)
    sizes = []
    matrix = np.zeros((l,l))
    filenames = list(files.keys())
    for i in range(l):
        for j in range(i+1, l):
            jacc,size = jaccard(files[filenames[i]], files[filenames[j]], k)
            matrix[i,j]  = jacc
            sizes.append(size)
            print(filenames[i], filenames[j], jacc)
    print(matrix)
    print(sizes)
    end = time.time()
    print('Temps écoulé : {}'.format(end-start))
    # Temps moyen : 260s à améliorer
    