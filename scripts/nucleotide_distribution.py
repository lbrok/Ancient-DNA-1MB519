import matplotlib.pyplot as plt
import pandas as pd

def nucleotide_distribution():
    # Reading fasta file to string
    fasta_file = open(r"C:/Users/46722/Documents/Tillämpad_bioinformatik/Ancient-DNA-1MB519/chromosomes/chr24.fa", "r")
    seq = fasta_file.read()
    fasta_file.close()

    # Reading file of intervals of ROIs
    intervals = pd.read_csv(r"C:/Users/46722/Documents/Tillämpad_bioinformatik/Ancient-DNA-1MB519/test_intervals.txt", sep = '\t')
    # Reading sequencing depth data
    data = pd.read_csv(
        r"C:\Users\46722\Documents\Tillämpad_bioinformatik\Ancient-DNA-1MB519\data\mammoth.chr24.average.depth.gz", sep='\t',
        comment='t', header=None,
        usecols=[1, 2, 3, 4])
    data.columns = ['pos', 'african_elephant', 'asian_elephant', 'woolly_mammoth']
    # Calculating nucleotide distributions in specified regions
    # print(len(seq))
    gc_counts = []
    depths = []

    for index, row in intervals.iterrows():
        start, stop = row['start'], row['stop']
        # print(f"Region interval: {start, stop}")
        region = seq[start:stop]
        no_ns = region.replace('N', '')    # Excluding positions with N
        gc_content = (region.count('G') + region.count('C'))/len(no_ns)
        gc_counts.append(gc_content)
        depths.append(data['woolly_mammoth'][start:stop].mean())

        # N_prop = region.count('N')/len(region)
        # N_counts.append(N_prop)
        # print(f"G: {G_prop} \nC: {C_prop} \nA: {A_prop} \nT: {T_prop} \nN: {N_prop}")

    plt.scatter(gc_counts, depths)
    plt.xlabel('GC-content')
    plt.ylabel('Sequence coverage depth')
    plt.show()

if __name__ == '__main__':
    nucleotide_distribution()