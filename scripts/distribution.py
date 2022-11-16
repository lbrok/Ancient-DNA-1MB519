import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import poisson

def plot_data():
    data = pd.read_csv(r"C:\Users\46722\Documents\Tillämpad_bioinformatik\Ancient-DNA-1MB519\data\mammoth.chr24.depth.gz", sep='\t', comment='t', header=None,
                       usecols=[1, 6, 7, 8, 9, 10])  # change in usecols for which columns you want to use
    data.columns = ['pos', 'mammoth_1', 'mammoth_2', 'mammoth_3', 'mammoth_4', 'mammoth_5']
    df = pd.DataFrame(data)
    df = df.iloc[3037418:]  # Remove telomeres
    subplot_nr=1
    fig = plt.figure()
    plt.rcParams["figure.figsize"] = (15, 12)
    fig.subplots_adjust(hspace=0.5, wspace=0.9)
    fig.suptitle("Histograms of sequencing coverage for all Woolly mammoths after preprocessing")
    for ind in ['mammoth_1', 'mammoth_2', 'mammoth_3', 'mammoth_4', 'mammoth_5']:
        mean = df[ind].mean()
        std = df[ind].std()
        print(mean, std)
        df2 = df[df[ind] < (mean+(2*std))]
        bins = max(df2[ind])
        plt.subplot(2, 3, subplot_nr)
        plt.hist(df2[ind], bins=bins)
        plt.title(f"{ind}", fontsize=10)
        plt.ylim(0, 2300000)
        plt.xlim(0, 70)
        plt.xlabel('sequencing coverage', fontsize=8)
        subplot_nr += 1
    plt.show()

if __name__ == '__main__':
    plot_data()