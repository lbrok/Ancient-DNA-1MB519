import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def plot_file():
    data = pd.read_csv('mammoth.chr24.depth.gz', sep='\t', comment='t', header=None,
                       usecols=[1, 6])  # change in usecols for which colums you want to use
    # data.columns = ['chromosome', 'chromosome position',
    #               'Asian-elephant', 'Asian-elephant',
    #              'African-elephant', 'African-elephant',
    #             'Woolly-mammoth-1',
    #            'Woolly-mammoth', 'Woolly-mammoth',
    #           'Woolly-mammoth', 'Woolly-mammoth']
    data.columns = ['pos', 'woolly_1']
    df = pd.DataFrame(data)

    # fig = plt.plot(df['pos'], df['woolly_1'], label='Woolly_Mammoth_1')
    # fig = plt.plot(df['pos'], df['woolly_2'], label='Woolly_Mammoth_2')
    # fig = plt.plot(df['pos'], df['woolly_3'], label='Woolly_Mammoth_3')
    # fig = plt.plot(df['pos'], df['woolly_4'], label='Woolly_Mammoth_4')
    # fig = plt.plot(df['pos'], df['woolly_5'], label='Woolly_Mammoth_5')

    plt.rcParams["figure.figsize"] = (12, 7)
    # plt.subplot(1, 2, 1)
    plt.hist(df['woolly_1'], bins=max(df['woolly_1']))
    plt.title('Histogram of sequencing coverage of Woolly mammoth 1')
    plt.xlabel('sequencing coverage')
    # plt.xlim(0, 100)
    plt.yscale('log')
    plt.show()

    # plt.subplot(1, 2, 2)
    # # df.drop(df.loc[df['woolly_1'] == 0].index, inplace=True)
    # df = df.iloc[2999995:]
    # plt.hist(df['woolly_1'], bins=max(df['woolly_1']))
    # plt.title('Histogram of sequencing coverage \nin interval [0:100] of Woolly mammoth 1\n not including telomere')
    # plt.xlabel('sequencing coverage')
    # # plt.xlim(0, 100)
    # plt.yscale('log')
    # plt.show()

if __name__ == '__main__':
    plot_file()
