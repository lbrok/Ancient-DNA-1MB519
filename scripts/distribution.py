import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import poisson

def fit_data():
    data = pd.read_csv('mammoth.chr24.depth.gz', sep='\t', comment='t', header=None,
                       usecols=[1, 6])  # change in usecols for which columns you want to use
    data.columns = ['pos', 'woolly_1']
    df = pd.DataFrame(data)
    df = df.iloc[3037500:]  # Remove telomeres
    mean = df['woolly_1'].mean()
    std = df['woolly_1'].std()
    print(mean, std)
    df2 = df[df['woolly_1'] < (mean+(2*std))]
    bins = max(df2['woolly_1'])
    plt.hist(df2['woolly_1'], bins=bins)
    plt.title('Histogram of sequencing coverage of Woolly mammoth 1, after preprocessing')
    plt.xlabel('sequencing coverage')
    plt.show()

if __name__ == '__main__':
    fit_data()