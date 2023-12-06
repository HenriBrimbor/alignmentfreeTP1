from loading import load_directory
from kmers import stream_kmers, kmer2str
import numpy as np
import time
import sys

code = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}


def jaccard(kmers_A, kmers_B):
    """
    fileA / fileB : lists of lists containing sequences on different chromosomes in A or B
    k : integer for the length of kmers
    """

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
    
    del kmers_A
    del kmers_B
    
    return intersect / union



if __name__ == "__main__":
    start = time.time()
    
    
    files = load_directory("data") # Load all the files in a dictionary
    k = 21                         # setting the length of kmers
    l = len(files)
    matrix       = np.zeros((l,l))
    filenames    = list(files.keys())
    kmers_files  = {}
    max_memory   = 0
    files_memory = np.sum([np.sum([sys.getsizeof(i) for i in file]) for file in files]) 

    # Rather than calculating every kmers of the same file at each iteration
    # we store them in a dictionnary outside of jaccard to calculate them only once
    for i in range(l):
        # if the kmers of file i are met for the first time build them
        if i not in kmers_files:
            kmers_files[i] = {}
            for list_Ai in files[filenames[i]]:
                for kmer, kmer_comp in stream_kmers(list_Ai, k):
                    kmers_files[i][kmer]      = kmers_files[i].get(kmer, 0) + 1
                    kmers_files[i][kmer_comp] = kmers_files[i].get(kmer_comp, 0) + 1
            
        for j in range(i+1, l):
            # if the kmers of file j are met for the first time build them
            if j not in kmers_files:
                    kmers_files[j] = {}
                    for list_Bi in files[filenames[j]]:
                        for kmer, kmer_comp in stream_kmers(list_Bi, k):
                            kmers_files[j][kmer]      = kmers_files[j].get(kmer, 0) + 1
                            kmers_files[j][kmer_comp] = kmers_files[j].get(kmer_comp, 0) + 1
                            
            matrix[i,j] = jaccard( kmers_files[i], kmers_files[j] )
            #print(filenames[i], filenames[j], matrix[i,j])

        max_memory = max(max_memory, np.sum( [sys.getsizeof(kmers_files[i]) for i in kmers_files.keys()]) )
        # All calculations done for file i, releasing the memory
        del kmers_files[i]
        
    print(matrix)
    end = time.time()
    
    # Amélioration : from 280s to 180s but uses more memory
    print('Time spent : {}'.format( end-start) )
    print('Max memory used with dict : {}'.format ( max_memory  ) )
    print('Max memory used with files : {}'.format( files_memory) )
    print('Sum : {} kb = {} Tb'.format( files_memory+max_memory, round((files_memory+max_memory)/10**9, 4) ) )
    
    with open('./README.md', 'w') as f:
        for index,i in enumerate(files.keys()):
            f.write( '#{} {}\n'.format(index, i))
        for line in matrix:
            string = ""
            for item in line:
                string += str(item) + '\t'
            string = string[:-2]+'\n'
            f.write(string)
        
    
    
    
    
    
    
    
    