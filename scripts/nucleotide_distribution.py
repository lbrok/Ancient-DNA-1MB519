# Reading fasta file to string
fasta_file = open(r"C:/Users/46722/Documents/Tillämpad_bioinformatik/Ancient-DNA-1MB519/chr24.fa", "r")
seq = fasta_file.read()
fasta_file.close()

import pandas as pd
# Reading text file of intervals of ROIs
# text_file = open(r"C:/Users/46722/Documents/Tillämpad_bioinformatik/Ancient-DNA-1MB519/test_intervals.txt", "r")
# intervals = text_file.read()
# text_file.close()

intervals = pd.read_csv(r"C:/Users/46722/Documents/Tillämpad_bioinformatik/Ancient-DNA-1MB519/test_intervals.txt", sep = '\t')

# Calculating nucleotide distributions in specified regions
#print(len(seq))
for index, row in intervals.iterrows():
    start, stop = row['start'], row['stop']
    print(f"Region interval: {start, stop}")
    region = seq[start:stop]
    no_Ns = region.replace('N', '')    # Excluding positions with N
    G_prop = region.count('G')/len(no_Ns)
    C_prop = region.count('C')/len(no_Ns)
    A_prop = region.count('A')/len(no_Ns)
    T_prop = region.count('T')/len(no_Ns)
    N_prop = region.count('N')/len(region)
    print(f"G: {G_prop} \nC: {C_prop} \nA: {A_prop} \nT: {T_prop} \nN: {N_prop}")
