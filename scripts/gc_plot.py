import matplotlib.pyplot as plt
import pandas as pd

def plot_gc():
    # Reading fasta file to string
    fasta_file = open(r"C:/Users/46722/Documents/Tillämpad_bioinformatik/Ancient-DNA-1MB519/chromosomes/chr24.fa", "r")
    seq = fasta_file.read()
    fasta_file.close()

    # Reading sequencing depth data
    data = pd.read_csv(
        r"C:\Users\46722\Documents\Tillämpad_bioinformatik\Ancient-DNA-1MB519\data\mammoth.chr24.depth.gz", sep='\t',
        comment='t', header=None,
        usecols=[1, 2, 3, 6])
    data.columns = ['pos', 'african_elephant', 'asian_elephant', 'woolly_mammoth']

    # Calculating depth and GC content in specified regions
    # print(len(seq))
    gc_counts = []
    depths = []
    telomere_length = 3037418
    threshold = 31.2    # Calculated based on genome average and std deviation

    for interval in list(range(telomere_length+1, len(seq), 1000)):     # Window size 1000 bp
        start, stop = interval, (interval + 1000)
        depth = data['woolly_mammoth'][start:stop].mean()
        #asian_depth = data['asian_elephant'][start:stop].mean()
        #african_depth = data['african_elephant'][start:stop].mean()
        asian_depth = 0
        african_depth = 0
        if depth <= threshold:
            if (asian_depth <= depth) or (african_depth <= depth):
                depths.append(depth)
                region = seq[start:stop]
                no_ns = region.replace('N', '')  # Excluding positions with N
                gc_content = (region.count('G') + region.count('C')) / len(no_ns)
                gc_counts.append(gc_content)

    # Plotting results
    plt.scatter(gc_counts, depths, s=5)
    plt.xlabel('GC-content')
    plt.ylabel('Sequence coverage depth')
    plt.show()

if __name__ == '__main__':
    plot_gc()