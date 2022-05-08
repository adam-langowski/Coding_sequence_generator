This program produces all possible coding sequences from a given amino acid sequence.
In addition, it calculates its moll mass, as well as GC content and hamming distance between all created DNA sequence. 

Example of needed amino_acids.txt file:

SPW

Output (codons.txt file):

Mass moll of given peptide sequence: 424.4483

All possible sequences based on given amino acids: 

1) AGTCCATGG
{'A': 2, 'C': 2, 'T': 2, 'G': 3}

GC content in this sequence: 0.56

Complement sequence: TCAGGTACC
...

{1: 'AGTCCATGG', 2: 'AGTCCGTGG', 3: 'AGTCCCTGG', 4: 'AGTCCTTGG', 5: 'AGCCCATGG', 6: 'AGCCCGTGG', 7: 'AGCCCCTGG', 8: 'AGCCCTTGG', 9: 'TCACCATGG', 10: 'TCACCGTGG', 11: 'TCACCCTGG', 12: 'TCACCTTGG', 13: 'TCCCCATGG', 14: 'TCCCCGTGG', 15: 'TCCCCCTGG', 16: 'TCCCCTTGG', 17: 'TCGCCATGG', 18: 'TCGCCGTGG', 19: 'TCGCCCTGG', 20: 'TCGCCTTGG', 21: 'TCTCCATGG', 22: 'TCTCCGTGG', 23: 'TCTCCCTGG', 24: 'TCTCCTTGG'}

