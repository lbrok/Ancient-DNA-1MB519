import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

path = '../Data/mammoth.chr24.depth'
start = 3037418
'''
incrementalStatistics(file_path, start)
A function to calculate the mean value and standard deviation while dividing the data and working in chunks.
RETURNS: A list, index 0 is the length of the file in the number of nucleotide, index 1 is the mean value, and index 2 is the std.
EXAMPLES: incrementalStatistics('../Data/mammoth.chr24.depth',3037418) = [24086192, 15.205945721930114, 14.45353032871932]
'''
def incrementalStatistics(filePath, start):
    n = 0 #Number of nucleotides
    m = 0 #Mean value
    s = 0 #Standard deviation
    for chunk in pd.read_csv(filePath, sep='\t', comment='t', header=None, usecols=[6,7,8,9,10], skiprows=start, chunksize=5000000):
        chunk = chunk.mean(axis=1) #Calculates the mean value of each input row
        for index, mammothMean in chunk.items(): #Iterates through each row
            n += 1 #Adds one nucleotide in the counter
            mOld = m #Stores the recent mean value
            m = m + ((mammothMean-m)/n) #Calculates an incremental mean value
            s = s + ((mammothMean-mOld)*(mammothMean-m)) #Calculates an incremental std
    
    s = (s/(n-1))**(0.5) #Calculates the final standard deviation
    return [n,m,s] #[Length, mean, std]

print(incrementalStatistics(path,start))