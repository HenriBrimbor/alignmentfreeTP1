
def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for i in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)


def stream_kmers(text, k):
    letters = ['A', 'C', 'T', 'G']
    mask = (1 << (2*k))-1 # décallage de 2 bits fois k nucléotides (pour k-mer codés en 2-bits)
    # pour k=2 encodage de 2 bits: 1<<4 -1 = 000001<<4 -1 = 010000 -1 = 001111
    # l'opérateur & ET logique va donc filtrer les 4 derniers bits uniquement donc les deux derniers nucléotides

    # initialisation des entiers, à voir comme des cases mémoire avec les bits
    kmer = 0
    comp = 0
    retenue = 0
    for index, letter in enumerate(text):
        if letter in letters:
            nucl = letters.index(letter)   # on cherche l'index de la lettre dans la liste, qui correspond à l'encodage 0,1,2,3 équivalent en binaire (0,0) (0,1) (1,0) (1,1)
            kmer <<= 2                     # on décale les bits vers la gauche par rapport au cadre de lecture
            kmer += nucl                   # on ajoute l'index composé de deux bits après le décalage donc sur la droite des précédents
            kmer &= mask                   # on applique le masque pour ne prendre que les k derniers caractères (prise en compte des k*2 derniers bits)

            # On fait la même chose avec le kmer complémentaire
            nucl_comp = (nucl+2) % 4
            if index+1-retenue < k:
                comp += nucl_comp<<(2*(index-retenue+1))
            # pour les k-1 premiers il ne faut rien renvoyer, car le masque prend les 0 à gauche en mémoire comme des A, alors qu'ils n'existent pas dans la séquence
            else: # index-retenue+1 >= k:
                comp += nucl_comp<<(2*(k))
                comp >>= 2
                yield kmer, comp
        else:
            # Dans le cas où l'on aurait une lettre qui n'est pas dans la liste des lettres, dans les k premières lettres
            # on doit l'ignorer mais décrémenter l'index
            if index+1 < k:
                retenue += 1

if __name__ == "__main__":
    k = 3
    for i in stream_kmers('AACACTG', k):
        print( kmer2str(i[0], k), kmer2str(i[1], k) )