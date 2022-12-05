import matplotlib.pyplot as plt
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os
import numpy as np
from window_slider import Slider
from Bio.Seq import Seq


# str(Seq('CGG').reverse_complement())


def plot_gc():
    outfile = str(input('Output file name: '))
    # Calculating depth and GC content in specified regions using these dictionaries
    # Dictionaries have keys from 0 to 100

    codons = ['CGA', 'CGG', 'CGT', 'CGC', 'GCA', 'GCG', 'GCT', 'GCC', 'TCT', 'TCC', 'TTA', 'TTG',
              'CTA', 'CTG', 'ACA', 'ACG', 'CTT', 'CTC', 'GTT', 'GTC', 'CCA', 'CCG', 'CCT', 'CCC',
              'GTA', 'GTG', 'TAA', 'TAG', 'TAT', 'AAT', 'AAC', 'GAA', 'GAG', 'TTT', 'TTC', 'TAC',
              'AAA', 'AAG', 'GAT', 'GAC', 'ATG', 'ATT', 'ACT', 'TGA', 'TGG', 'TGT', 'TGC', 'ACC',
              'ATA', 'CAA', 'CAG', 'CAT', 'CAC', 'AGA', 'AGG', 'AGT', 'AGC', 'TCA', 'TCG', 'GGA',
              'GGG', 'GGT', 'GGC', 'ATC']

    print(len(codons))

    non_codons = ['NNN', 'NNG', 'NNC', 'NNA', 'NNT', 'NGG', 'NGC', 'NGA', 'NGT', 'NCC', 'NCG', 'NCA', 'NCT',
                  'NAA', 'NAG', 'NAC', 'NAT', 'NTT', 'NTG', 'NTC', 'NTA', 'GNN', 'CNN', 'ANN', 'TNN', 'GGN',
                  'GCN', 'GAN', 'GTN', 'CCN', 'CGN', 'CAN', 'CTN', 'AAN', 'AGN', 'ACN', 'ATN', 'TTN', 'TGN',
                  'TCN', 'TAN']
    print(len(non_codons))
    # codons_dict = {i: [0, 0, 0] for i in codons}
    codons_dict = {k: {i: 0 for i in np.around(np.arange(-1, 1, 0.001), 2)} for k in codons}
    # mOld_dict = {i: 0 for i in codons}

    # Defining variables
    window_size = 1000  # Set window size for sliding window analysis
    mammoth_average = 13.64328801191864  # Pooled mammoth genome average
    elephant_average = 33.158036502652195  # Pooled elephant genome average
    chunk_size = 1000000  # Size of data batches being read, in rows at a time
    names = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
             'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25',
             'chr26', 'chr27', 'chrX']
    telomeres = [3000000, 3000000, 0, 3000000, 3000000, 3000000, 0, 3000000, 3000000, 3000000, 3000000, 3000000,
                 3000000,
                 3000000,
                 300000, 3000000, 3000000, 3000000, 3000000, 3000000, 0, 3000000,
                 3000000, 3000000, 3000000, 3000000, 3000000, 3000000]

    bucket_size = 3
    overlap_count = 1

    with open(
            r"C:\Users\miche\PycharmProjects\datachromo\Data\LoxAfr4_DQ188829.fa") as handle:
        for values in SimpleFastaParser(
                handle):  # Loads the Fasta header and the corresponding sequence in a list [header,sequence]
            if values[
                0] in names:  # If the fasta header is present in the names list it should be included in the # calculations
                chr_nr = names.index(values[0])  # Getting the index of the fasta header in the names list
                for idx, data in enumerate(pd.read_csv(
                        r"C:\Users\miche\PycharmProjects\datachromo\Data\mammoth." + names[
                            chr_nr] + ".depth", sep='\t',
                        comment='t', header=None,
                        usecols=[2, 3, 4, 5, 6, 7, 8, 9], chunksize=chunk_size,
                        skiprows=telomeres[chr_nr])):
                    # Excluding the 5th mammoth
                    data.columns = ['asian_elephant_1', 'asian_elephant_2', 'african_elephant_1', 'african_elephant_2',
                                    'woolly_mammoth_1', 'woolly_mammoth_2', 'woolly_mammoth_3', 'woolly_mammoth_4']

                    # Computing average depth at each position of each species, removing old columns from dataframe
                    data['average_elephant'] = (data['asian_elephant_1'] + data['asian_elephant_2'] + data[
                        'african_elephant_1']
                                                + data['african_elephant_2']) / 4
                    data.drop(['asian_elephant_1', 'asian_elephant_2', 'african_elephant_1', 'african_elephant_2'],
                              axis='columns', inplace=True)
                    data['average_mammoth'] = (data['woolly_mammoth_1'] + data['woolly_mammoth_2'] + data[
                        'woolly_mammoth_3']
                                               + data['woolly_mammoth_4']) / 4
                    data.drop(['woolly_mammoth_1', 'woolly_mammoth_2', 'woolly_mammoth_3', 'woolly_mammoth_4'],
                              axis='columns', inplace=True)

                    for interval in list(range(0, len(data), window_size)):  # Window size 1000 bp
                        start, stop = interval, (interval + window_size)

                        mammoth_depth = data['average_mammoth'].iloc[start:stop].mean()
                        elephant_depth = data['average_elephant'].iloc[start:stop].mean()
                        slider = Slider(bucket_size, overlap_count)

                        if mammoth_depth <= 28.7812709234281:  # Removing outliers (value from (average + 2*std)) 13.655611525936575+2*7.5628296987457775
                            if mammoth_depth != 0 or elephant_depth != 0:
                                ratio = ((mammoth_depth / mammoth_average) - (elephant_depth / elephant_average)) / \
                                        ((mammoth_depth / mammoth_average) + (elephant_depth / elephant_average))
                                ratio = round(ratio, 2)

                                slider.fit(
                                    np.array(list(values[1][telomeres[chr_nr] + idx * chunk_size + start:telomeres[
                                                                                                             chr_nr] + idx * chunk_size + stop])))
                                new_slide = True
                                if ratio > 1:
                                    print(f'Positions: {[start, stop]} Ratio: {ratio}]')
                                while new_slide:
                                    window_data = slider.slide()
                                    window_data = ''.join(window_data)
                                    if window_data in non_codons or len(window_data) < 3:
                                        break
                                    # Updating the dictionary for the correct codon, by adding 1 occurrance
                                    codons_dict[window_data][ratio] = (codons_dict[window_data][ratio]) + 1

                        # print(
                        #  f'Chromosome {names[chr_nr]} {round(((telomeres[chr_nr] + idx * chunk_size + start) / len(values[1])) * 100, 2)}%')
                data = 0
                print(f'Chromosome {names[chr_nr]} DONE')

            values = 0

        # ratio = []
        # occurrence = []
        # Save results to file
        #print(codons_dict)
        print('Writing and plotting')
        for i in codons:
            f = open(outfile + i + ".txt", "w")
            ratio_counts = []
            codon_occurrence = []
            for j in codons_dict[i]:
                if codons_dict[i][j] != 0:
                    ratio_counts.append(j)
                    codon_occurrence.append(codons_dict[i][j])
                f.write(str(float(j)) + '\t' + str(codons_dict[i][j]) + '\n')
            # ratio.append(ratio_counts)
            # print(ratio)
            # occurrence.append(codon_occurrence)
            # print(occurence)
            f.close()

        # # Fit trend line to data
        # # Parameters from the fit of the polynomial of the codon
        # p1 = np.polyfit(occurrence[0], ratio[0], deg=1)
        #
        # # Fit trend line to data
        # # Parameters from the fit of the polynomial of the reverse-complement codon
        # p2 = np.polyfit(occurrence[1], ratio[1], deg=1)
        #
        # # Model the data using the parameters of the fitted straight line of the codon
        # y_model_1 = np.polyval(p1, occurrence[0])
        #
        # # Model the data using the parameters of the fitted straight line of the reverse-complement codon
        # y_model_2 = np.polyval(p2, occurrence[1])
        #
        # # Mean codon
        # y_bar_1 = np.mean(ratio[0])
        # # Mean reverse codon
        # y_bar_2 = np.mean(ratio[1])
        # # Coefficient of determination, R²
        # R2_1 = np.sum((y_model_1 - y_bar_1) ** 2) / np.sum((ratio[0] - y_bar_1) ** 2)
        # # Coefficient of determination, R²
        # R2_2 = np.sum((y_model_2 - y_bar_2) ** 2) / np.sum((ratio_counts[1] - y_bar_2) ** 2)
        #
        # # Plotting results
        # plt.scatter(occurrence[0], ratio[0], s=5, color='blue')
        # plt.axhline(0, color='r', linewidth=0.5)
        # plt.xlabel('Occurences of "CGG" codon in window')
        # plt.title(f'Window size {window_size} bp')
        # plt.ylabel('Mammoth depth ratio - Elephant depth ratio / \nMammoth depth ratio + Elephant depth ratio')
        #
        # # Line of best fit
        # xlim = plt.xlim()
        # plt.plot(np.array(xlim), p1[1] + p1[0] * np.array(xlim), label=f'Line of Best Fit, R² = {R2_1:.2f}', color='black')
        #
        # # Plotting results
        # plt.scatter(occurrence[1], ratio[1], s=5, color='orange')
        #
        #
        # # Line of best fit
        # xlim = plt.xlim()
        # plt.plot(np.array(xlim), p2[1] + p2[0] * np.array(xlim), label=f'Line of Best Fit, R² = {R2_2:.2f}',
        #          color='black')
        # plt.legend(fontsize=8)
        # plt.show()


if __name__ == '__main__':
    plot_gc()
