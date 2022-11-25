import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def plot_gc():
    # Reading fasta file to string
    fasta_file = open(r"C:/Users/46722/Documents/Applied_bioinformatics/Ancient-DNA-1MB519/chromosomes/chr24.fa", "r")
    seq = fasta_file.read()
    fasta_file.close()

    # Reading sequencing depth data
    data = pd.read_csv(
        r"C:\Users\46722\Documents\Applied_bioinformatics\Ancient-DNA-1MB519\data\mammoth.chr24.depth.gz", sep='\t',
        comment='t', header=None,
        usecols=[1, 2, 3, 4, 5, 6, 7, 8, 9])    # Excluding the 5th mammoth
    data.columns = ['pos', 'asian_elephant_1', 'asian_elephant_2', 'african_elephant_1', 'african_elephant_2',
                    'woolly_mammoth_1', 'woolly_mammoth_2', 'woolly_mammoth_3', 'woolly_mammoth_4']

    # Computing average depth at each position of each species, removing old columns from dataframe
    data['average_elephant'] = (data['asian_elephant_1'] + data['asian_elephant_2'] + data['african_elephant_1']
                                + data['african_elephant_2']) / 4
    data.drop(['asian_elephant_1', 'asian_elephant_2', 'african_elephant_1', 'african_elephant_2'],
              axis='columns', inplace=True)
    data['average_mammoth'] = (data['woolly_mammoth_1'] + data['woolly_mammoth_2'] + data['woolly_mammoth_3']
                               + data['woolly_mammoth_4']) / 4
    data.drop(['woolly_mammoth_1', 'woolly_mammoth_2', 'woolly_mammoth_3', 'woolly_mammoth_4'],
              axis='columns', inplace=True)

    # Calculating depth and GC content in specified regions
    gc_dict = {i: [] for i in range(100)}
    # g_dict = {i: [] for i in range(100)}
    # c_dict = {i: [] for i in range(100)}
    ratios = []

    # Defining variables
    window_size = 1000  # Set window size for sliding window analysis
    telomere_length = 3037418   # Known value, true for all except chromosomes 3, 7 and 21 (no telomere sequence)
    mammoth_average = 13.34943     # Pooled mammoth genome average
    elephant_average = 30   # Pooled elephant genome average

    for interval in list(range(telomere_length+1, len(seq), window_size)):     # Window size 1000 bp
        start, stop = interval, (interval + window_size)
        mammoth_depth = data['average_mammoth'][start:stop].mean()
        elephant_depth = data['average_elephant'][start:stop].mean()

        if mammoth_depth <= 28.7812709:  # Removing outliers (value from (average + 2*std))
            if mammoth_depth == 0 or elephant_depth == 0:
                continue
            else:
                ratio = ((mammoth_depth/mammoth_average)-(elephant_depth/elephant_average))/\
                        ((mammoth_depth/mammoth_average)+(elephant_depth/elephant_average))
                ratios.append(ratio)
                region = seq[start:stop]
                no_ns = region.replace('N', '')  # Excluding positions with N
                gc_content = ((region.count('G') + region.count('C')) / len(no_ns))*100   # Calculating GC-content
                gc_dict[round(gc_content)].append(ratio)
                # g_dict[round((region.count('G')/len(no_ns))*100)].append(ratio)
                # c_dict[round((region.count('C')/len(no_ns))*100)].append(ratio)

    # Calculating average depth in each GC-window
    ratio_counts = []
    for i in gc_dict:
        if gc_dict[i] != []:
            ratio_counts.append(sum(gc_dict[i]) / len(gc_dict[i]))
        else:
            ratio_counts.append(None)

    # Calculating average depth in each G-window
    # ratio_counts_g = []
    # for i in g_dict:
    #     if g_dict[i] != []:
    #         ratio_counts_g.append(sum(g_dict[i]) / len(g_dict[i]))
    #     else:
    #         ratio_counts_g.append(None)
    #
    # # Calculating average depth in each C-window
    # ratio_counts_c = []
    # for i in c_dict:
    #     if c_dict[i] != []:
    #         ratio_counts_c.append(sum(c_dict[i]) / len(c_dict[i]))
    #     else:
    #         ratio_counts_c.append(None)

    # Fit trend line to data

    # Parameters from the fit of the polynomial
    ratio_counts_no_none = [i for i in ratio_counts if i != None]
    x = [i for i in range(len(ratio_counts_no_none))]
    p = np.polyfit(x, ratio_counts_no_none, deg=1)

    # Model the data using the parameters of the fitted straight line
    y_model = np.polyval(p, x)

    # Mean
    y_bar = np.mean(ratio_counts_no_none)
    # Coefficient of determination, R²
    R2 = np.sum((y_model - y_bar) ** 2) / np.sum((ratio_counts_no_none - y_bar) ** 2)

    print(f'R² = {R2:.2f}')

    # Plotting results
   #  plt.subplot(3, 1, 1)
    plt.scatter([i for i in range(100)], ratio_counts, s=5)
    # plt.axhline(0, color='r', linewidth=0.5)
    plt.xlabel('GC-content')
    plt.ylabel('Mammoth depth ratio/Elephant depth ratio')
    # plt.subplot(3, 1, 2)
    # plt.scatter([i for i in range(100)], ratio_counts_g, s=5)
    # plt.axhline(0, color='r', linewidth=0.5)
    # plt.xlabel('G-content')
    # plt.subplot(3, 1, 3)
    # plt.scatter([i for i in range(100)], ratio_counts_c, s=5)
    # plt.axhline(0, color='r', linewidth=0.5)
    # plt.xlabel('C-content')

    # Line of best fit
    xlim = plt.xlim()
    plt.plot(np.array(xlim), p[1] + p[0] * np.array(xlim), label=f'Line of Best Fit, R² = {R2:.2f}', color='black')
    plt.legend(fontsize=8)
    plt.show()


if __name__ == '__main__':
    plot_gc()
