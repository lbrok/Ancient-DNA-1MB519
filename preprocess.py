import pandas as pd
import csv
from tabulate import tabulate


def preprocess(inputfile, telo_start, telo_end, chr, std, mean):
    result = []
    data = pd.read_csv(inputfile, sep='\t', comment='t', header=None,
                       usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], skiprows=[i for i in range(telo_start, telo_end)],
                       chunksize=100000)
    for chunk in data:
        chunk.columns = ['chr', 'Pos', 'Asian-1', 'Asian-2', 'African-1', 'African-2', 'Woolly-mammoth-1',
                         'Woolly-mammoth-2', 'Woolly-mammoth-3', 'Woolly-mammoth-4', 'Woolly-mammoth-5']
        df = pd.DataFrame(chunk)

        df['avr_pooled_asian'] = (df['Asian-1'] + df['Asian-2']) / 2
        df.drop(['Asian-1', 'Asian-2'], axis='columns', inplace=True)
        df['avr_pooled_african'] = (df['African-1'] + df['African-2']) / 2
        df.drop(['African-1', 'African-2'], axis='columns', inplace=True)

        df['avr_pooled_mammoth'] = (df['Woolly-mammoth-1'] + df['Woolly-mammoth-2'] + df['Woolly-mammoth-3']
                                    + df['Woolly-mammoth-4'] + df['Woolly-mammoth-5']) / 5

        df.drop(['Woolly-mammoth-1', 'Woolly-mammoth-2', 'Woolly-mammoth-3', 'Woolly-mammoth-4', 'Woolly-mammoth-5'],
                axis='columns', inplace=True)

        # std = df[['avr_pooled_mammoth']].std()

        outliers = mean + (2 * std)
        # print(outliers)

        df2 = df[df['avr_pooled_mammoth'] < outliers]
        # print(df2)
        result.append(df2)
        # print(result)

    result = pd.concat(result)

    content = tabulate(result.values.tolist(), tablefmt="plain")

    with open(f"{chr}_preprocessed.bed", "w") as outfile:
        outfile.write(content)


if __name__ == '__main__':
    preprocess('mammoth.chr24.depth.gz', 0, 3000000, 'chr24', 14, 7)
