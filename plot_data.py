import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def plot_file():
    data = pd.read_csv('mammoth.chr24.depth', sep='\t', comment='t', header=None,
                       usecols=[1, 6, 7, 8, 9, 10])  # change in usecols for which colums you want to use
    # data.columns = ['chromosome', 'chromosome position',
    #               'Asian-elephant', 'Asian-elephant',
    #              'African-elephant', 'African-elephant',
    #             'Woolly-mammoth-1',
    #            'Woolly-mammoth', 'Woolly-mammoth',
    #           'Woolly-mammoth', 'Woolly-mammoth']
    data.columns = ['pos', 'wolly_1', 'wolly_2', 'wolly_3', 'wolly_4', 'wolly_5']
    df = pd.DataFrame(data)

    fig = plt.plot(df['pos'], df['wolly_1'], label='Woolly_Mammoth_1')
    fig = plt.plot(df['pos'], df['wolly_2'], label='Woolly_Mammoth_2')
    fig = plt.plot(df['pos'], df['wolly_3'], label='Woolly_Mammoth_3')
    fig = plt.plot(df['pos'], df['wolly_4'], label='Woolly_Mammoth_4')
    fig = plt.plot(df['pos'], df['wolly_5'], label='Woolly_Mammoth_5')
    plt.title('showing coverage corresponding to the position')
    plt.legend()
    plt.xlabel('position')
    plt.ylabel('coverage')
    plt.show()


if __name__ == '__main__':
    plot_file()
