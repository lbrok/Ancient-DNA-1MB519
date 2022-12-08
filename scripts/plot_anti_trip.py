import matplotlib.pyplot as plt
import os
import pandas as pd
from Bio.Seq import Seq
import math


def plot_gc():
    codons = [('CGA', str(Seq('CGA').reverse_complement())), ('CGG', str(Seq('CGG').reverse_complement())),
              ('CGT', str(Seq('CGT').reverse_complement())), ('CGC', str(Seq('CGC').reverse_complement())),
              ('GCA', str(Seq('GCA').reverse_complement())), ('GCT', str(Seq('GCT').reverse_complement())),
              ('GCC', str(Seq('GCC').reverse_complement())), ('TCT', str(Seq('TCT').reverse_complement())),
              ('TCC', str(Seq('TCC').reverse_complement())), ('TTA', str(Seq('TTA').reverse_complement())),
              ('TTG', str(Seq('TTG').reverse_complement())), ('CTA', str(Seq('CTA').reverse_complement())),
              ('CTG', str(Seq('CTG').reverse_complement())), ('ACA', str(Seq('ACA').reverse_complement())),
              ('CTT', str(Seq('CTT').reverse_complement())), ('CTC', str(Seq('CTC').reverse_complement())),
              ('GTT', str(Seq('GTT').reverse_complement())), ('GTC', str(Seq('GTC').reverse_complement())),
              ('CCA', str(Seq('CCA').reverse_complement())), ('CCT', str(Seq('CCT').reverse_complement())),
              ('CCC', str(Seq('CCC').reverse_complement())), ('GTA', str(Seq('GTA').reverse_complement())),
              ('GTG', str(Seq('GTG').reverse_complement())), ('TAT', str(Seq('TAT').reverse_complement())),
              ('AAT', str(Seq('AAT').reverse_complement())), ('GAA', str(Seq('GAA').reverse_complement())),
              ('TTT', str(Seq('TTT').reverse_complement())), ('GAT', str(Seq('GAT').reverse_complement())),
              ('ATG', str(Seq('ATG').reverse_complement())), ('ACT', str(Seq('ACT').reverse_complement())),
              ('TGA', str(Seq('TGA').reverse_complement())), ('ACC', str(Seq('GGT').reverse_complement()))]

    average_ratios = []
    stdev_ratios = []
    occurrences = []

    for i in range(len(codons)):

        trip = pd.read_csv(r"C:\Users\miche\PycharmProjects\ancient dna\output_" + codons[i][0] + ".txt",
                           sep='\t', comment='t', header=None)
        anti_trip = pd.read_csv(r"C:\Users\miche\PycharmProjects\ancient dna\output_" + codons[i][1] + ".txt",
                                sep='\t', comment='t', header=None)

        codon_occurrence = trip[1] + anti_trip[1]
        ratio_counts = trip[0]

        average_ratio = 0
        stdev_ratio = 0
        for i in range(len(trip[0])):
            average_ratio = average_ratio + (ratio_counts[i] * codon_occurrence[i])
            stdev_ratio = stdev_ratio + ((ratio_counts[i] - (average_ratio / sum(codon_occurrence))) ** 2)
        average_ratios.append(average_ratio / sum(codon_occurrence))
        stdev_ratios.append(math.sqrt(stdev_ratio / sum(codon_occurrence)))
        occurrences.append(sum(codon_occurrence))

        # Plotting results
    plot_color = '#1f7864'
    plt.subplot(2, 1, 1)
    plt.bar([str(f"{i}/{j}") for i, j in codons], occurrences, color=plot_color)
    plt.xticks(rotation=90)
    plt.title(f"Window size 1000 bp")
    plt.ylabel(f"Number of occurrences")
    plt.subplot(2, 1, 2)
    plt.errorbar([str(f"{i}/{j}") for i, j in codons], average_ratios, stdev_ratios, linestyle='None', marker='o', color=plot_color)
    plt.xticks(rotation=90)
    plt.grid(axis='x', color='whitesmoke')
    plt.ylabel(f"Average depth ratio")
    plt.show()


if __name__ == '__main__':
    plot_gc()
