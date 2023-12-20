from loading import load_directory
from kmers import stream_kmers, kmer2str
import numpy as np
import time
import sys
import heapq
from math import inf

code = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

def triFusion(L):
    if len(L) == 1:
        return L
    else:
        return fusion( triFusion(L[:len(L)//2]) , triFusion(L[len(L)//2:]) )

def fusion(A,B):
    if len(A) == 0:
        return B
    elif len(B) == 0:
        return A
    elif A[0] <= B[0]:
        return [A[0]] + fusion( A[1:] , B )
    else:
        return [B[0]] + fusion( A , B[1:] )

def XorShift(x):
    x ^= (x << 13) & 0xFFFFFFFFFFFFFFFF
    x ^= (x >> 7)  & 0xFFFFFFFFFFFFFFFF
    x ^= (x << 17) & 0xFFFFFFFFFFFFFFFF
    return x

def sampling( sequence, k, s ):
    sketch = [+inf] * s
    for kmer in stream_kmers(sequence, k):
        idx = sketch.index(max(sketch))
        if kmer < sketch[idx]:
            sketch[idx] = kmer
    return sketch

def better_sampling( file, k, s ):
    sketch = [-inf] * s
    heapq.heapify( sketch )
    for sequence in file:
        for kmer in stream_kmers(sequence, k):
            # on prend la plus petite valeur du sketch == la plus grande avec inversement de signe
            elem = heapq.nsmallest(1, sketch)

            # si l'opposé du kmer choisi est plus grand alors on remplace l'élément par le kmer
            # == si le kmer choisi est plus petit que le kmer max
            # ie : plus grand que le kmer min  // opposé // plus petit que le kmer max
            if -kmer > elem:
                heapq.heappop(sketch)      # on enlève la plus petite valeur == la plus grande avec inversion de signe
                heapq.push(sketch, -kmer)
    return sketch

def better_sampling_Xorshift( fna_file, k, s ):
    sketch = [-inf] * s
    heapq.heapify( sketch )
    for sequence in fna_file:
        for kmer in stream_kmers(sequence, k):
            xor = XorShift(kmer)
            
            # si le kmer choisi est plus petit alors on remplace l'élément par le kmer (xorshifté)
            # on compare avec la plus petite valeur du sketch == la plus grande avec inversement de signe
            if -xor > heapq.nsmallest(1, sketch)[0]:
                #heapq.heappop(sketch)  # on enlève la plus petite valeur == la plus grande avec inversion de signe
                #heapq.heappush(sketch, -xor)  # on insère le -xor plus grand == xor plus petit
                heapq.heapreplace(sketch, -xor)
    return sketch

def jaccard_list(file_A, file_B):
        """
        Function that returns the jaccard value of two sorted lists of kmers A and B by using iterative comparaison
        """
        i, j   = 0, 0
        l1, l2 = len(file_A), len(file_B)
        Union  = 0
        Inter  = 0
        
        while i!=l1 and j!=l2:
            Union += 1
            if file_A[i] == file_B[j]:
                Inter += 1
                i += 1
                j += 1
            else:
                if file_A[i] < file_B[j]:
                    i += 1
                else:
                    j += 1

        return Inter / Union


if __name__ == "__main__":
    start = time.time()
    
    try:
        k = int(sys.argv[1])
    except:
        raise ValueError('k not defined')
    try:
        s = int(sys.argv[2])
    except:
        raise ValueError('s not defined')
    print('k:', k)
    print('s:', s)
    
    files = load_directory("data") # Load all the files in a dictionary
    
    l            = len(files)
    matrix       = np.zeros((l,l))
    filenames    = list(files.keys())
    files_sizes  = np.sum([np.sum([sys.getsizeof(i) for i in file]) for file in files])
    files_memory = 0 
    kmers_files  = {}
    # Rather than calculating every kmers of the same file at each iteration
    # we store them in a dictionnary outside of jaccard to calculate them only once
    for i in range(len(files)):
        for j in range(i+1, l):
            if i not in kmers_files:
                kmers_files[i] = np.sort( better_sampling_Xorshift( files[filenames[i]], k, s) )
            if j not in kmers_files:
                kmers_files[j] = np.sort( better_sampling_Xorshift( files[filenames[j]], k, s) )
            matrix[i,j] = jaccard_list( kmers_files[i], kmers_files[j] )
        
        files_memory = max(files_memory, sys.getsizeof(kmers_files) )
    print(matrix)
    end = time.time()
    

    print('Time spent : {}'.format( end-start) )
    print('Max memory used with files : {}'.format ( files_sizes   ) )
    print('Max memory used with dict  : {}'.format( files_memory ) )
    print('Sum : {} b = {} Tb'.format( files_memory+files_sizes, round((files_memory+files_sizes)/10**6, 5) ) )
    
    with open('./README.md', 'w') as f:
        f.write('Parameters used :\n')
        f.write('s : {} \n'.format( s ) )
        f.write('k : {} \n'.format( k ) )
        f.write('time run : {}s\n'.format(end-start))
        f.write('memory used : {} b = {} Mb\n'.format( files_memory+files_sizes, round((files_memory+files_sizes)/10**6, 5) ) )
        for index,i in enumerate(files.keys()):
            f.write( '#{} {}\n'.format(index, i))
        for line in matrix:
            string = ""
            for item in line:
                string += str(item) + '\t'
            string = string[:-2]+'\n'
            f.write(string)
        f.write('\n'*2)
        f.write('Obtaining approximated results of the previous ones\nHowever the heapq method seems to exponentially extend the process time')
        
    
    
    
    
    
    
    
    