import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

telomeres = [3037418,3037418,0,3037418,3037418,3037418,0,3037418,3037418,3037418,3037418,3037418,3037418,3037418,3037418,3037418,3037418,3037418,3037418,3037418,0,3037418,3037418,3037418,3037418,3037418,3037418,3037418]
#names = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chrX']
names = ['chr24']
path = 'D:/Data/mammoth.'

'''
incrementalStatistics(file_path, start)
A function to calculate the mean value and standard deviation while dividing the data and working in chunks.
RETURNS: A list, index 0 is the length of the file in the number of nucleotide, index 1 is the mean value, and index 2 is the std.
EXAMPLES: incrementalStatistics('../Data/mammoth.chr24.depth',3037418) = [24086192, 15.205945721930114, 14.45353032871932]
'''
def incrementalStatistics(filePath, names, telomeres):
    n = 0 #Number of nucleotides
    m = 0 #Mean value
    s = 0 #Standard deviation
    #meanValue = 33.108845127365235 #Used when removing outliers, this is the previous average
    #stdValue = 10.777236719119456 #Used when removing outliers, this is the previous average
    for i in range(len(names)):
        for chunk in pd.read_csv(filePath+names[i]+'.depth', sep='\t', comment='t', header=None, usecols=[2,3,4,5], skiprows=telomeres[i], chunksize=10000000):
            chunk = chunk.mean(axis=1) #Calculates the mean value of each input row
            for index, mammothMean in chunk.items(): #Iterates through each row
                #if mammothMean < meanValue + 2*stdValue: #Add / remove if statement if outliers should be removed or not
                n += 1 #Adds one nucleotide in the counter
                mOld = m #Stores the recent mean value
                m = m + ((mammothMean-m)/n) #Calculates an incremental mean value
                s = s + ((mammothMean-mOld)*(mammothMean-m)) #Calculates an incremental std
        print(f'Chromosome {i+1} of {len(telomeres)}')
    s = (s/(n-1))**(0.5) #Calculates the final standard deviation
    return [n,m,s] #[Length, mean, std]

print(incrementalStatistics(path,names,telomeres))