import matplotlib.pyplot as plt
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
import os
import numpy as np
from window_slider import Slider

def plot_gc():

        codon_occurrence = 

        # Fit trend line to data
        # Parameters from the fit of the polynomial
        p = np.polyfit(codon_occurrence, ratio_counts, deg=1)

        # Model the data using the parameters of the fitted straight line
        y_model = np.polyval(p, codon_occurrence)

        # Mean
        y_bar = np.mean(ratio_counts)
        # Coefficient of determination, R²
        R2 = np.sum((y_model - y_bar) ** 2) / np.sum((ratio_counts - y_bar) ** 2)

        # Plotting results
        plt.scatter(codon_occurrence, ratio_counts, s=5)
        plt.axhline(0, color='r', linewidth=0.5)
        plt.xlabel('Occurences of "AAA" codon in window')
        plt.title(f'Window size {window_size} bp')
        plt.ylabel('Mammoth depth ratio - Elephant depth ratio / \nMammoth depth ratio + Elephant depth ratio')

        # Line of best fit
        xlim = plt.xlim()
        plt.plot(np.array(xlim), p[1] + p[0] * np.array(xlim), label=f'Line of Best Fit, R² = {R2:.2f}', color='black')
        plt.legend(fontsize=8)
        plt.errorbar(codon_occurrence, ratio_counts, fmt='o')
        # yerr=std_counts
        plt.show()

if __name__ == '__main__':
    plot_gc()
