import matplotlib.pyplot as plt
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os
import numpy as np


def plot_gc():
    outfile=str(input('Output file name: '))
    # Calculating depth and GC content in specified regions using these dictionaries
    # Dictionaries have keys from 0 to 100, values are a list, index 0: avergae, 1: standard deviation, 2: number of datapoints
    gc_dict = {i: [0,0,0] for i in range(101)}
    #g_dict = {i: [0,0] for i in range(101)}
    #c_dict = {i: [0,0] for i in range(101)}
    mOld_dict = {i: 0 for i in range(101)}

    # Defining variables
    window_size = 10000  # Set window size for sliding window analysis
    mammoth_average = 13.34943063462263     # Pooled mammoth genome average
    elephant_average = 33.001500559085684   # Pooled elephant genome average
    chunk_size = 1200000 #Size of data batches being read, in rows at a time
    telomeres = [3000000,3000000,0,3000000,3000000,3000000,0,3000000,3000000,3000000,3000000,3000000,3000000,3000000,300000,3000000,3000000,3000000,3000000,3000000,0,3000000,3000000,3000000,3000000,3000000,3000000,3000000]
    names = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chrX']

    with open("D:/Data/LoxAfr4_DQ188829.fa", encoding="ascii") as handle:
        for values in SimpleFastaParser(handle): #Loads the Fasta header and the corresponding sequence in a list [header,sequence]
            if values[0] in names: #If the fasta header is present in the names list it should be included in the calculations
                chr_nr = names.index(values[0]) #Getting the index of the fasta header in the names list
                for idx,data in enumerate(pd.read_csv(
                    r"D:/Data/mammoth."+names[chr_nr]+".depth", sep='\t',
                    comment='t', header=None,
                    usecols=[2, 3, 4, 5, 6, 7, 8, 9], chunksize=chunk_size, skiprows=telomeres[chr_nr])):   # Excluding the 5th mammoth
                    data.columns = [ 'asian_elephant_1', 'asian_elephant_2', 'african_elephant_1', 'african_elephant_2',
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
                
                    for interval in list(range(0, len(data), window_size)):     # Window size 1000 bp
                        start, stop = interval, (interval + window_size)
                        mammoth_depth = data['average_mammoth'].iloc[start:stop].mean()
                        elephant_depth = data['average_elephant'].iloc[start:stop].mean()

                        if mammoth_depth <= 28.7812709234281:  # Removing outliers (value from (average + 2*std)) 13.655611525936575+2*7.5628296987457775
                            if mammoth_depth != 0 or elephant_depth != 0:
                                ratio = ((mammoth_depth/mammoth_average)-(elephant_depth/elephant_average))/\
                                    ((mammoth_depth/mammoth_average)+(elephant_depth/elephant_average))
                                region = values[1][telomeres[chr_nr]+idx*chunk_size+start:telomeres[chr_nr]+idx*chunk_size+stop]
                                no_ns = region.replace('N', '')  # Excluding positions with N
                                if len(no_ns) == 0:
                                    continue
                                else:
                                    gc_content = ((region.count('G') + region.count('C')) / len(no_ns))*100   # Calculating GC-content
                                    gc_dict_key = round(gc_content)

                                    #Saving the old average used for std calculations
                                    mOld_dict[gc_dict_key] = gc_dict[gc_dict_key][0]
                                    #Updating the averages of different GC%, and the number of datapoints
                                    #New average being the old average + (ratio-old average) / number of datapoints, including the new one
                                    gc_dict[gc_dict_key] = [gc_dict[gc_dict_key][0]+((ratio-gc_dict[gc_dict_key][0])/(gc_dict[gc_dict_key][2]+1)),gc_dict[gc_dict_key][1],(gc_dict[gc_dict_key][2])+1]
                
                                    #Updating the std (variance?)
                                    gc_dict[gc_dict_key][1] = gc_dict[gc_dict_key][1] + ((ratio-mOld_dict[gc_dict_key])*(ratio-gc_dict[gc_dict_key][0]))

                    print(f'Chromosome {names[chr_nr]} {round(((telomeres[chr_nr]+idx*chunk_size+start)/len(values[1]))*100,2)}%')
                data = 0
                print(f'Chromosome {names[chr_nr]} DONE')
        values = 0

    print('Writing and plotting')
    os.chdir("./GC_data")
    f = open(outfile+".txt", "a")
    # Calculating average depth in each GC-window
    GC_counts = []
    ratio_counts = []
    std_counts = []
    for i in gc_dict:
        if gc_dict[i] != [0,0,0] and gc_dict[i][2] != 1:
            GC_counts.append(i)
            gc_dict[i][1] = (gc_dict[i][1]/(gc_dict[i][2]-1))**(0.5) #Final std calc
            ratio_counts.append(gc_dict[i][0])
            std_counts.append(gc_dict[i][1]*2) #2stds for a 95% confidence interval
        f.write(str(round(float(gc_dict[i][0]),3))+'\t'+str(round(float(gc_dict[i][1]),3))+'\t'+str(int(gc_dict[i][2]))+'\n')
    f.close()

    # Fit trend line to data
    # Parameters from the fit of the polynomial
    p = np.polyfit(GC_counts,ratio_counts, deg=1)

    # Model the data using the parameters of the fitted straight line
    y_model = np.polyval(p, GC_counts)

    # Mean
    y_bar = np.mean(ratio_counts)
    # Coefficient of determination, R??
    R2 = np.sum((y_model - y_bar) ** 2) / np.sum((ratio_counts - y_bar) ** 2)

    # Plotting results
    plt.scatter(GC_counts, ratio_counts, s=5)
    plt.axhline(0, color='r', linewidth=0.5)
    plt.xlabel('GC-content [%]')
    plt.title(f'Window size {window_size} bp')
    plt.ylabel('Mammoth depth ratio - Elephant depth ratio / \nMammoth depth ratio + Elephant depth ratio')

    # Line of best fit
    xlim = plt.xlim()
    plt.plot(np.array(xlim), p[1] + p[0] * np.array(xlim), label=f'Line of Best Fit, R?? = {R2:.2f}', color='black')
    plt.legend(fontsize=8)
    plt.errorbar(GC_counts, ratio_counts, yerr=std_counts, fmt='o')
    plt.show()

if __name__ == '__main__':
    plot_gc()
