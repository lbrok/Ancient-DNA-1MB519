import matplotlib.pyplot as plt
import pandas as pd
import math

codons = ['CGA', 'CGG', 'CGT', 'CGC', 'GCA', 'GCG', 'GCT', 'GCC', 'TCT', 'TCC', 'TTA', 'TTG',
          'CTA', 'CTG', 'ACA', 'ACG', 'CTT', 'CTC', 'GTT', 'GTC', 'CCA', 'CCG', 'CCT', 'CCC',
          'GTA', 'GTG', 'TAA', 'TAG', 'TAT', 'AAT', 'AAC', 'GAA', 'GAG', 'TTT', 'TTC', 'TAC',
          'AAA', 'AAG', 'GAT', 'GAC', 'ATG', 'ATT', 'ACT', 'TGA', 'TGG', 'TGT', 'TGC', 'ACC',
          'ATA', 'CAA', 'CAG', 'CAT', 'CAC', 'AGA', 'AGG', 'AGT', 'AGC', 'TCA', 'TCG', 'GGA',
          'GGG', 'GGT', 'GGC', 'ATC']
codons.sort()


def average_plot():
    average_ratios = []
    stdev_ratios = []
    occurrences = []
    for codon in codons:
        data = pd.read_csv(r"C:\Users\46722\Documents\Applied_bioinformatics\Ancient-DNA-1MB519\codon_data\genome_wide" + codon + ".txt",
                                sep='\t', comment='t', header=None)
        codon_occurrence = data[1]
        ratio_counts = data[0]
        average_ratio = 0
        stdev_ratio = 0
        for i in range(len(data[0])):
            average_ratio = average_ratio + (ratio_counts[i]*codon_occurrence[i])
            stdev_ratio = stdev_ratio + ((ratio_counts[i] - (average_ratio/sum(data[1])))**2)
        average_ratios.append(average_ratio/sum(data[1]))
        stdev_ratios.append(math.sqrt(stdev_ratio/sum(data[1])))
        occurrences.append(sum(data[1]))

    # Plotting results
    plt.subplot(2, 1, 1)
    plt.bar(codons, occurrences, color='#1f7864')
    plt.xticks(rotation=90)
    plt.title(f"Window size 1000 bp", pad=20)
    plt.ylabel(f"Number of occurrences of trinucleotide genome wide", fontsize=8, labelpad=20)
    plt.subplot(2, 1, 2)
    plt.scatter(codons, average_ratios, linestyle='None', s=10, marker='o', color='#1f7864')
    plt.errorbar(codons, average_ratios, stdev_ratios, linestyle='None', fmt='', color='#1f7864')
    plt.xticks(rotation=90)
    plt.grid(axis='x', color='whitesmoke')
    plt.ylabel(f"Average depth ratio values for trinucleotide genome wide", fontsize=8, labelpad=5)
    plt.show()


if __name__ == '__main__':
    average_plot()
