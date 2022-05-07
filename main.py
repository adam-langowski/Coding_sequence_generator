# Adam Langowski 
# Temat projektu:
# Generator sekwencji kodujących (wej. sekw. aminokwasów,wyj. zbiór sekwencji kodujących)

import itertools

gencode_dic = {
    'A': ['GCA', 'GCC', 'CGC', 'GCT'],
    'C': ['TGC', 'TGT'],
    'D': ['GAC', 'GAT'],
    'E': ['GAA', 'GAG'],
    'F': ['TTC', 'TTT'],
    'G': ['GGA', 'GGC', 'GGG', 'GGT'],
    'H': ['CAC', 'CAT'],
    'I': ['ATA', 'ATT', 'ATC'],
    'K': ['AAA', 'AAG'],
    'L': ['CTA', 'CTG', 'CTC', 'CTT', 'TTA', 'TTG'],
    'M': ['ATG'],
    'N': ['AAC', 'AAT'],
    'P': ['CCA', 'CCG', 'CCC', 'CCT'],
    'Q': ['CCA', 'CAG'],
    'R': ['AAG', 'AGA', 'CGA', 'CGC', 'CGG', 'CGT'],
    'S': ['AGT', 'AGC', 'TCA', 'TCC', 'TCG', 'TCT'],
    'T': ['ACA', 'ACC', 'ACG', 'ACT'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'W': ['TGG'],
    'Y': ['TAC', 'TAT'],
    '*': ['TAA', 'TAG', 'TGA']
}

protein_weights = {
    'A': 89.0932,
    'C': 121.1582,
    'D': 133.1027,
    'E': 147.1293,
    'F': 165.1891,
    'G': 75.0666,
    'H': 155.1546,
    'I': 131.1729,
    'K': 146.1876,
    'L': 131.1729,
    'M': 149.2113,
    'N': 132.1179,
    'O': 255.3134,
    'P': 115.1305,
    'Q': 146.1445,
    'R': 174.201,
    'S': 105.0926,
    'T': 119.1192,
    'U': 168.0532,
    'V': 117.1463,
    'W': 204.2252,
    'Y': 181.1885
}


# Generator producing all possible codons from given amino acids sequence


def retranlation(pep):
    codons = [gencode_dic.get(amino_acid, 'XXX') for amino_acid in pep]
    for cod in itertools.product(*codons):
        yield "".join(cod)


# Function calculating moll mass of given peptide sequence


def moll_mass(pep):
    kDa = 0
    for char in pep:
        kDa += protein_weights.get(char, 0)
    print(f'Mass moll of given peptide sequence: {kDa:.4f}\n')


# Class with sequence created based on retranslation


class Sequence:
    def __init__(self, id, seq):
        self.id = id
        self.seq = seq

    def nucleotides(self, seq):  # Function counting amount of each nucleotide in retranslated sequence
        a = c = t = g = 0
        for char in seq:
            if char == 'A':
                a += 1
            elif char == 'C':
                c += 1
            elif char == 'T':
                t += 1
            elif char == 'G':
                g += 1
        d = {'A': a, 'C': c, 'T': t, 'G': g}
        return d

    def gc_content(self, seq):  # Function calculating GC content in created DNA sequence
        c_count = seq.count('C')
        g_count = seq.count('G')
        return f'GC content in this sequence: {float((c_count + g_count) / len(seq)):.2f}'

    def complement(self, seq):        # Function creating complement sequence
        d_complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        complement_seq = ''
        for char in seq:
            complement_seq += d_complement[char]
        return f'Complement sequence: {complement_seq}\n'


# Reading from file:

file = open('amino_acids.txt')
peptides = file.readline()
file.close()

# Output:

print(f'Your input sequence: {peptides}\n')
moll_mass(peptides)  # Displaying moll mass

print(f'All possible sequences based on given amino acids: \n')

seq_list = {}
count = 1
for cod in retranlation(peptides):
    s = Sequence(id=count, seq=cod)  # creating codon object
    seq_list[s.id] = s.seq
    print(f'{count}) {cod}')
    count += 1
    print(s.nucleotides(s.seq))
    print(s.gc_content(s.seq))
    print(s.complement(s.seq))


def hamming_distance(seq1, seq2):  # Function calculating hamming distances between all sequences in created dic
    return sum([1 for n1, n2 in zip(seq1, seq2) if n1 != n2])


distances = []
distance = 0
for seq1, seq2 in itertools.combinations(seq_list.keys(), 2):
    distance = hamming_distance(seq_list.get(seq1), seq_list.get(seq2))
    distances.append(distance)
    print(f'Hamming distance beetwen {seq1} and {seq2}: {distance}')
print(f'\nMinimum hamming distance between two sequences: {min(distances)}')

# Writing to file
file2 = open('codons.txt', 'w')
file2.write(str(seq_list))
file2.close()