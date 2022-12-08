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
    mOld_dict = {i: 0 for i in range(101)}
    remove_chunks = 0

    # Defining variables
    window_size = 100  # Set window size for sliding window analysis
    neanderthal_average = 13.336119102787015    # Pooled mammoth genome average
    sapien_average = 33.10835482336245  # Pooled elephant genome average
    chunk_size = 5000000 #Size of data batches being read, in rows at a time
    
    telomeres = [10000,10000,60000,10000,10000,60000,10000,10000,10000,60000,60000,60000,19020000,19000000,20000000,60000,0,10000,60000,60000,9411193,16050000]
    fasta_names = ['1 dna:chromosome chromosome:GRCh37:1:1:249250621:1', '2 dna:chromosome chromosome:GRCh37:2:1:243199373:1',
    '3 dna:chromosome chromosome:GRCh37:3:1:198022430:1', '4 dna:chromosome chromosome:GRCh37:4:1:191154276:1',
    '5 dna:chromosome chromosome:GRCh37:5:1:180915260:1', '6 dna:chromosome chromosome:GRCh37:6:1:171115067:1',
    '7 dna:chromosome chromosome:GRCh37:7:1:159138663:1', '8 dna:chromosome chromosome:GRCh37:8:1:146364022:1',
    '9 dna:chromosome chromosome:GRCh37:9:1:141213431:1', '10 dna:chromosome chromosome:GRCh37:10:1:135534747:1',
    '11 dna:chromosome chromosome:GRCh37:11:1:135006516:1', '12 dna:chromosome chromosome:GRCh37:12:1:133851895:1',
    '13 dna:chromosome chromosome:GRCh37:13:1:115169878:1', '14 dna:chromosome chromosome:GRCh37:14:1:107349540:1',
    '15 dna:chromosome chromosome:GRCh37:15:1:102531392:1', '16 dna:chromosome chromosome:GRCh37:16:1:90354753:1',
    '17 dna:chromosome chromosome:GRCh37:17:1:81195210:1', '18 dna:chromosome chromosome:GRCh37:18:1:78077248:1',
    '19 dna:chromosome chromosome:GRCh37:19:1:59128983:1', '20 dna:chromosome chromosome:GRCh37:20:1:63025520:1',
    '21 dna:chromosome chromosome:GRCh37:21:1:48129895:1', '22 dna:chromosome chromosome:GRCh37:22:1:51304566:1']
    depth_names = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

    with open("D:/Data/Neanderthal/homo-sapiens-reference.fa", encoding="ascii") as handle:
        for values in SimpleFastaParser(handle): #Loads the Fasta header and the corresponding sequence in a list [header,sequence]
            if values[0] in fasta_names: #If the fasta header is present in the names list it should be included in the calculations
                chr_nr = fasta_names.index(values[0]) #Getting the index of the fasta header in the names list
                for i in range(int(telomeres[chr_nr]/chunk_size)): #Loop that calculates how many chunks that can be skipped. This saves memory instead of using skiprows in read csv
                    telomeres[chr_nr] -= chunk_size
                    remove_chunks += 1
                for idx,data in enumerate(pd.read_csv(
                    r"D:/Data/Neanderthal/neanderthal-human."+depth_names[chr_nr]+".depth.gz", sep='\t',
                    comment='t', header=None,
                    usecols=[2, 3, 4, 5], chunksize=chunk_size, skiprows=telomeres[chr_nr])):   # Excluding the 5th mammoth
                    data.columns = ['neanderthal_1', 'neanderthal_2', 'sapien_1', 'sapien_2']
                    if remove_chunks == 0:
                        # Computing average depth at each position of each species, removing old columns from dataframe
                        data['average_neanderthal'] = (data['neanderthal_1'] + data['neanderthal_2']) / 2
                        data.drop(['neanderthal_1', 'neanderthal_2'], axis='columns', inplace=True)
                        data['average_sapien'] = (data['sapien_1'] + data['sapien_2']) / 2
                        data.drop(['sapien_1', 'sapien_2'],axis='columns', inplace=True)
                    
                        for interval in list(range(0, len(data), window_size)):     # Window size 1000 bp
                            start, stop = interval, (interval + window_size)
                            neanderthal_depth = data['average_neanderthal'].iloc[start:stop].mean()
                            sapien_depth = data['average_sapien'].iloc[start:stop].mean()

                            if neanderthal_depth <= 90.743432296371:  # Removing outliers (value from (average + 2*std))
                                if neanderthal_depth != 0 or sapien_depth != 0:
                                    ratio = ((neanderthal_depth/neanderthal_average)-(sapien_depth/sapien_average))/\
                                        ((neanderthal_depth/neanderthal_average)+(sapien_depth/sapien_average))
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
                        print(f'Chromosome {depth_names[chr_nr]} {round(((telomeres[chr_nr]+idx*chunk_size+start)/len(values[1]))*100,2)}%')
                    else:
                        remove_chunks -= 1
                data = 0
                print(f'Chromosome {depth_names[chr_nr]} DONE')
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

if __name__ == '__main__':
    plot_gc()
