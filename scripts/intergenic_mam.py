import matplotlib.pyplot as plt
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os
import numpy as np


def plot_gc():
    outfile = str(input('Output file name: '))
    # Calculating depth and GC content in specified regions using these dictionaries
    # Dictionaries have keys from 0 to 100, values are a list, index 0: avergae, 1: standard deviation, 2: number of datapoints
    # genic_regions = []
    # non_genic_regions = []

    genic_gc_dict = {i: [0, 0, 0] for i in range(101)}
    genic_mOld_dict = {i: 0 for i in range(101)}

    non_gc_dict = {i: [0, 0, 0] for i in range(101)}
    non_mold_dict = {i: 0 for i in range(101)}

    # Defining variables
    window_size = 10000  # Set window size for sliding window analysis
    mammoth_average = 13.34943063462263  # Pooled mammoth genome average
    elephant_average = 33.001500559085684  # Pooled elephant genome average
    chunk_size = 1000000  # Size of data batches being read, in rows at a time
    telomeres = [3000000, 3000000, 0, 3000000, 3000000, 3000000, 0, 3000000, 3000000, 3000000, 3000000, 3000000,
                 3000000, 3000000, 300000, 3000000, 3000000, 3000000, 3000000, 3000000, 0, 3000000, 3000000, 3000000,
                 3000000, 3000000, 3000000, 3000000]
    telomeres_1 = [3000000, 3000000, 0, 3000000, 3000000, 3000000, 0, 3000000, 3000000, 3000000, 3000000, 3000000,
                   3000000, 3000000, 300000, 3000000, 3000000, 3000000, 3000000, 3000000, 0, 3000000, 3000000, 3000000,
                   3000000, 3000000, 3000000, 3000000]

    # 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
    # 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25',
    # 'chr26', 'chr27', 'chrX'
    names = ['chr24']

    with open(
            r"C:\Users\miche\PycharmProjects\datachromo\Data\LoxAfr4_DQ188829.fa") as handle:
        for values in SimpleFastaParser(
                handle):  # Loads the Fasta header and the corresponding sequence in a list [header,sequence]
            if values[
                0] in names:  # If the fasta header is present in the names list it should be included in the calculations
                chr_nr = names.index(values[0])  # Getting the index of the fasta header in the names list
                chromosome = names[chr_nr]
                print(chromosome)

                genes = pd.read_csv(
                    r"C:\Users\miche\PycharmProjects\datachromo\Data\Loxodonta-Africana_reference_genes.gff.txt",
                    comment='"', sep='\s+', header=None, usecols=[0, 3, 4])
                genes.columns = ['chr', 'start', 'stop']
                genes = genes[genes['chr'] == chromosome]
                genes = genes.reset_index()
                non_genes = pd.DataFrame(columns=['chr', 'start', 'stop'])
                old_row = None
                for index, row in genes.iterrows():
                    if old_row is None:
                        new_row = pd.Series(data={'chr': names[chr_nr], 'start': 0, 'stop': row[2] - 1}, name=index)

                    else:
                        new_row = pd.Series(data={'chr': names[chr_nr], 'start': old_row[3] + 1, 'stop': row[2] - 1},
                                            name=index)
                    non_genes = non_genes.append(new_row, ignore_index=False)
                    old_row = row
                print(genes)
                number = 0

                # Loop that calculates how many chunks that can be skipped. This saves memory instead of using skiprows in read csv
                for index, region in genes.iterrows():
                    # print('skiprows:', region[2] - telomeres[chr_nr])

                    data = pd.read_csv(
                        r"C:\Users\miche\PycharmProjects\datachromo\Data\mammoth." + names[
                            chr_nr] + ".depth", sep='\t',
                        comment='t', header=None,
                        usecols=[2, 3, 4, 5, 6, 7, 8, 9], nrows=(region[3] - region[2]),
                        skiprows=(region[2]))  # Excluding the 5th mammoth

                    print('skiprows:', region[2])
                    print('data:', region[3] - region[2])

                    data.columns = ['asian_elephant_1', 'asian_elephant_2', 'african_elephant_1',
                                    'african_elephant_2',
                                    'woolly_mammoth_1', 'woolly_mammoth_2', 'woolly_mammoth_3', 'woolly_mammoth_4']

                    # Computing average depth at each position of each species, removing old columns from dataframe
                    data['average_elephant'] = (data['asian_elephant_1'] + data['asian_elephant_2'] + data[
                        'african_elephant_1']
                                                + data['african_elephant_2']) / 4
                    data.drop(
                        ['asian_elephant_1', 'asian_elephant_2', 'african_elephant_1', 'african_elephant_2'],
                        axis='columns', inplace=True)
                    data['average_mammoth'] = (data['woolly_mammoth_1'] + data['woolly_mammoth_2'] + data[
                        'woolly_mammoth_3']
                                               + data['woolly_mammoth_4']) / 4
                    data.drop(['woolly_mammoth_1', 'woolly_mammoth_2', 'woolly_mammoth_3', 'woolly_mammoth_4'],
                              axis='columns', inplace=True)

                    # # Looping through genic regions

                    mammoth_depth = data['average_mammoth'].mean()
                    elephant_depth = data['average_elephant'].mean()
                    number += 1
                    print(number, mammoth_depth, elephant_depth)

                    if mammoth_depth <= 28.7812709234281:  # Removing outliers (value from (average + 2*std)) 13.655611525936575+2*7.5628296987457775
                        if mammoth_depth != 0 or elephant_depth != 0:
                            ratio = ((mammoth_depth / mammoth_average) - (elephant_depth / elephant_average)) / \
                                    ((mammoth_depth / mammoth_average) + (elephant_depth / elephant_average))
                            print(ratio)
                            if ratio > 1:
                                print(f'Positions: {[start, stop]} Ratio: {ratio}]')
                            region = values[1][region[2]:region[3]]

                            no_ns = region.replace('N', '')  # Excluding positions with N
                            if len(no_ns) == 0:
                                continue
                            else:

                                gc_content = ((region.count('G') + region.count('C')) / len(
                                    no_ns)) * 100  # Calculating GC-content
                                genic_gc_dict_key = round(gc_content)

                                # Saving the old average used for std calculations
                                genic_mOld_dict[genic_gc_dict_key] = genic_gc_dict[genic_gc_dict_key][0]
                                # Updating the averages of different GC%, and the number of datapoints
                                # New average being the old average + (ratio-old average) / number of datapoints, including the new one
                                genic_gc_dict[genic_gc_dict_key] = [genic_gc_dict[genic_gc_dict_key][0] + (
                                        (ratio - genic_gc_dict[genic_gc_dict_key][0]) / (
                                        genic_gc_dict[genic_gc_dict_key][2] + 1)),
                                                                    genic_gc_dict[genic_gc_dict_key][1],
                                                                    (genic_gc_dict[genic_gc_dict_key][2]) + 1]

                                # Updating the std (variance?)
                                genic_gc_dict[genic_gc_dict_key][1] = genic_gc_dict[genic_gc_dict_key][1] + (
                                        (ratio - genic_mOld_dict[genic_gc_dict_key]) * (
                                        ratio - genic_gc_dict[genic_gc_dict_key][0]))

                # Looping through intergenic regions

                # for index, region in non_genes.iterrows():  # Window size 1000 bp
                #   start, stop = region[1], region[2]           # fel med detta!!!
                #  print(start, stop)
                # mammoth_depth = data['average_mammoth'].iloc[start:stop].mean()
                # elephant_depth = data['average_elephant'].iloc[start:stop].mean()
                # print(mammoth_depth)

                # if mammoth_depth <= 28.7812709234281:  # Removing outliers (value from (average + 2*std)) 13.655611525936575+2*7.5628296987457775
                #     if mammoth_depth != 0 or elephant_depth != 0:
                #         ratio = ((mammoth_depth / mammoth_average) - (elephant_depth / elephant_average)) / \
                #                 ((mammoth_depth / mammoth_average) + (elephant_depth / elephant_average))
                #         if ratio > 1:
                #             print(f'Positions: {[start, stop]} Ratio: {ratio}]')
                #         region = values[1][telomeres[chr_nr] + idx * chunk_size + start:telomeres[chr_nr] + idx * chunk_size + stop]
                #
                #         no_ns = region.replace('N', '')  # Excluding positions with N
                #         if len(no_ns) == 0:
                #             continue
                #         else:
                #             gc_content = ((region.count('G') + region.count('C')) / len(
                #                 no_ns)) * 100  # Calculating GC-content
                #             non_gc_dict_key = round(gc_content)
                #
                #             # Saving the old average used for std calculations
                #             non_mold_dict[non_gc_dict_key] = non_gc_dict[non_gc_dict_key][0]
                #             # Updating the averages of different GC%, and the number of datapoints
                #             # New average being the old average + (ratio-old average) / number of datapoints, including the new one
                #             non_gc_dict[non_gc_dict_key] = [non_gc_dict[non_gc_dict_key][0] + (
                #                     (ratio - non_gc_dict[non_gc_dict_key][0]) / (
                #                         non_gc_dict[non_gc_dict_key][2] + 1)),
                #                                                 non_gc_dict[non_gc_dict_key][1],
                #                                                 (non_gc_dict[non_gc_dict_key][2]) + 1]
                #
                #             # Updating the std (variance?)
                #             non_gc_dict[non_gc_dict_key][1] = non_gc_dict[non_gc_dict_key][1] + (
                #                     (ratio - non_mold_dict[non_gc_dict_key]) * (
                #                         ratio - non_gc_dict[non_gc_dict_key][0]))

        data = 0
        # print(f'Chromosome {names[chr_nr]} DONE')
    values = 0

    print('Writing and plotting')
    f = open(outfile + '_genic' + ".txt", "w")
    # Calculating average depth in each GC-window

    for i in genic_gc_dict:
        if genic_gc_dict[i] != [0, 0, 0] and genic_gc_dict[i][2] != 1:
            genic_gc_dict[i][1] = (genic_gc_dict[i][1] / (genic_gc_dict[i][2] - 1)) ** (0.5)  # Final std calc
        f.write(
            str(round(float(genic_gc_dict[i][0]), 3)) + '\t' + str(round(float(genic_gc_dict[i][1]), 3)) + '\t' + str(
                int(genic_gc_dict[i][2])) + '\n')
    f.close()

    print('Writing and plotting')
    f = open(outfile + '_non_genic' + ".txt", "w")
    # Calculating average depth in each GC-window
    #
    # for i in non_gc_dict:
    #     if non_gc_dict[i] != [0, 0, 0] and non_gc_dict[i][2] != 1:
    #         non_gc_dict[i][1] = (non_gc_dict[i][1] / (non_gc_dict[i][2] - 1)) ** (0.5)  # Final std calc
    #     f.write(
    #         str(round(float(non_gc_dict[i][0]), 3)) + '\t' + str(round(float(non_gc_dict[i][1]), 3)) + '\t' + str(
    #             int(non_gc_dict[i][2])) + '\n')
    # f.close()


if __name__ == '__main__':
    plot_gc()
