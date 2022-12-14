import matplotlib.pyplot as plt
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
import copy
import os
import numpy as np


def plot_gc():
    outfile = str(input('Output file name: '))
    genic_regions = []
    non_genic_regions = []

    # Defining variables
    window_size = 10000  # Set window size for sliding window analysis
    mammoth_average = 13.34943063462263  # Pooled mammoth genome average
    elephant_average = 33.001500559085684  # Pooled elephant genome average
    chunk_size = 30000000  # Size of data batches being read, in rows at a time
    # telomeres = [3000000, 3000000, 0, 3000000, 3000000, 3000000, 0, 3000000, 3000000, 3000000, 3000000, 3000000,
    #              3000000, 3000000, 3000000, 3000000, 3000000, 3000000, 3000000, 3000000, 0, 3000000, 3000000, 3000000,
    #              3000000, 3000000, 3000000, 3000000]
    # names = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
    #          'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25',
    #          'chr26', 'chr27', 'chrX']
    telomeres = [3000000]
    names = ['chr1']

    with open(
            r"C:\Users\46722\Documents\Applied_bioinformatics\Ancient-DNA-1MB519\data\Loxodonta-Africana_reference.fa") as handle:
        for values in SimpleFastaParser(handle):  # Loads the Fasta header and the corresponding sequence in a list [header,sequence]
            if values[0] in names:  # If the fasta header is present in the names list it should be included in the calculations
                chr_nr = names.index(values[0])  # Getting the index of the fasta header in the names list
                chromosome = names[chr_nr]

                # Creating dataframes for gene intervals
                genes = pd.read_csv(
                    r"C:\Users\46722\Documents\Applied_bioinformatics\Ancient-DNA-1MB519\data\Loxodonta-Africana_reference_genes.gff.txt",
                    comment='"', sep='\t', header=None, usecols=[0, 3, 4])
                genes.columns = ['chr', 'start', 'stop']
                genes = genes[genes['chr'] == chromosome]
                genes = genes.reset_index()

                # Creating dataframe for non-gene intervals
                non_genes = pd.DataFrame(columns=['chr', 'start', 'stop'])
                old_row = None
                for index, row in genes.iterrows():
                    if old_row is None:
                        new_row = pd.Series(
                            data={'chr': chromosome, 'start': genes['stop'].values[0] + 1, 'stop': len(values[1])},
                            name=index)

                    else:
                        new_row = pd.Series(data={'chr': chromosome, 'start': row[3] + 1, 'stop': old_row[2] - 1},
                                            name=index)
                    non_genes = non_genes.append(new_row, ignore_index=False)
                    old_row = copy.deepcopy(row)
                non_genes = non_genes.append(
                    pd.Series(data={'chr': chromosome, 'start': 0, 'stop': old_row[2] - 1}, name=len(values[1])),
                    ignore_index=True)

                for idx, data in enumerate(pd.read_csv(
                        r"C:\Users\46722\Documents\Applied_bioinformatics\Ancient-DNA-1MB519\data/mammoth." + names[
                            chr_nr] + ".depth.gz", sep='\t',
                        comment='t', header=None,
                        usecols=[2, 3, 4, 5, 6, 7, 8, 9], chunksize=chunk_size,
                        skiprows=telomeres[chr_nr])):  # Excluding the 5th mammoth
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

                    # Looping through genic regions
                    for index, region in genes.iterrows():
                        start, stop = (region[2] - telomeres[chr_nr]), (region[3] - telomeres[chr_nr])
                        mammoth_depth = data['average_mammoth'].iloc[start:stop].mean()
                        elephant_depth = data['average_elephant'].iloc[start:stop].mean()

                        # if mammoth_depth <= 28.7812709234281:  # Removing outliers (value from (average + 2*std)) 13.655611525936575+2*7.5628296987457775
                        if mammoth_depth != 0 or elephant_depth != 0:
                            ratio = ((mammoth_depth / mammoth_average) - (elephant_depth / elephant_average)) / \
                                    ((mammoth_depth / mammoth_average) + (elephant_depth / elephant_average))
                            if ratio > 1:
                                print(f'Positions: {[start, stop]} Ratio: {ratio}]')
                            region = values[1][telomeres[chr_nr] + idx * chunk_size + start:telomeres[chr_nr] + idx * chunk_size + stop]
                            no_ns = region.replace('N', '')  # Excluding positions with N
                            if len(no_ns) == 0:
                                continue
                            else:
                                genic_regions.append(ratio)

                    #Looping through intergenic regions
                    for index, region in non_genes.iterrows():  # Window size 1000 bp
                        start, stop = region[1], region[2]
                        mammoth_depth = data['average_mammoth'].iloc[start:stop].mean()
                        elephant_depth = data['average_elephant'].iloc[start:stop].mean()

                        if mammoth_depth <= 28.7812709234281:  # Removing outliers (value from (average + 2*std)) 13.655611525936575+2*7.5628296987457775
                            if mammoth_depth != 0 or elephant_depth != 0:
                                ratio = ((mammoth_depth / mammoth_average) - (elephant_depth / elephant_average)) / \
                                        ((mammoth_depth / mammoth_average) + (elephant_depth / elephant_average))
                                if ratio > 1:
                                    print(f'Positions: {[start, stop]} Ratio: {ratio}]')
                                region = values[1][telomeres[chr_nr] + idx * chunk_size + start:telomeres[chr_nr] + idx * chunk_size + stop]
                                no_ns = region.replace('N', '')  # Excluding positions with N
                                if len(no_ns) == 0:
                                    continue
                                else:
                                    non_genic_regions.append(ratio)

                data = 0
                print(f'Chromosome {names[chr_nr]} DONE')
        values = 0

    # Save results to file
    print('Writing and plotting')
    f = open(outfile + '_genic' + ".txt", "a")
    print(len(genic_regions))
    f.write(str(genic_regions))
    f.close()
    e = open(outfile + '_non_genic' + ".txt", "a")
    print(len(non_genic_regions))
    e.write(str(non_genic_regions))
    e.close()

    # Plotting results
    # plt.scatter(genic_regions, col='blue', s=5)
    # plt.scatter(non_genic_regions, col='blue', s=5)
    # plt.xlabel('GC-content [%]')
    # plt.title(f'Window size {window_size} bp')
    # plt.ylabel('Mammoth depth ratio - Elephant depth ratio / \nMammoth depth ratio + Elephant depth ratio')
    # plt.subplot(3, 1, 3)
    # plt.scatter([i for i in range(101)], ratio_counts_c, s=5)
    # plt.axhline(0, color='r', linewidth=0.5)
    # plt.xlabel('C-content [%]')

    plt.show()


if __name__ == '__main__':
    plot_gc()
