import matplotlib.pyplot as plt
import os
import pandas as pd

codons = ['CGA', 'CGG', 'CGT', 'CGC', 'GCA', 'GCG', 'GCT', 'GCC', 'TCT', 'TCC', 'TTA', 'TTG',
          'CTA', 'CTG', 'ACA', 'ACG', 'CTT', 'CTC', 'GTT', 'GTC', 'CCA', 'CCG', 'CCT', 'CCC',
          'GTA', 'GTG', 'TAA', 'TAG', 'TAT', 'AAT', 'AAC', 'GAA', 'GAG', 'TTT', 'TTC', 'TAC',
          'AAA', 'AAG', 'GAT', 'GAC', 'ATG', 'ATT', 'ACT', 'TGA', 'TGG', 'TGT', 'TGC', 'ACC',
          'ATA', 'CAA', 'CAG', 'CAT', 'CAC', 'AGA', 'AGG', 'AGT', 'AGC', 'TCA', 'TCG', 'GGA',
          'GGG', 'GGT', 'GGC', 'ATC']


def plot_gc():
        for codon in codons:
                data = pd.read_csv(r"C:\Users\46722\Documents\Applied_bioinformatics\Ancient-DNA-1MB519\codon_data/output" + codon + ".txt",
                                        sep='\t', comment='t', header=None)
                codon_occurrence = data[1]
                ratio_counts = data[0]

                # Plotting results
                plt.scatter(ratio_counts, codon_occurrence, s=5)
                plt.axhline(0, color='r', linewidth=0.5)
                plt.xlabel('Mammoth depth ratio - Elephant depth ratio / \nMammoth depth ratio + Elephant depth ratio')
                plt.title(f'Window size 1000 bp')
                plt.yscale('log')
                plt.ylabel(f"Occurrences of codons genome wide")

                # Save plot to figure
                os.chdir('../codon_data')
                plt.savefig(f"genome_wide_{codon}.png")


if __name__ == '__main__':
    plot_gc()