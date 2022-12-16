import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def readFiles(file_path):
    #file_path = input('File path to data: ')
    data = pd.read_csv(file_path, sep='\t', comment='t', header=None)
    data.columns = ['Mean', 'Std', 'NumOfWindows']
    data.insert(0, 'GC%', [i for i in range(101)])
    return data

def meanCalc(data):
    sum_NumOfWindows = data.sum()['NumOfWindows']
    #Mean value calculations
    mean_GC = 0
    for i in range(len(data)):
        mean_GC += i * data['NumOfWindows'][i]
    mean_GC = mean_GC/sum_NumOfWindows
    return mean_GC

def stdCalc(data, mean_GC):
    sum_NumOfWindows = data.sum()['NumOfWindows']
    #Standard deviation calculations
    variance = 0
    for i in range(len(data)):
        variance += ((i-mean_GC)**2) * data['NumOfWindows'][i]
    variance = variance/sum_NumOfWindows
    std_GC = variance**0.5
    return std_GC

def findOutliers(data, mean, std, minNumWindows=False):
    data['Outlier'] = ""
    for i in range(len(data)):
        if minNumWindows != False:
            if data.loc[i,'NumOfWindows'] < minNumWindows:
                data.loc[i,'Outlier'] = True
            else:
                data.loc[i,'Outlier'] = False
        else:
                if i < round(mean-numOfStds*std) or i > round(mean+numOfStds*std):
                    data.loc[i,'Outlier'] = True
                else:
                    data.loc[i,'Outlier'] = False

def dataFrameToList(data, column):
    lst = data[column].where(data['Outlier']==False)
    lst = list(lst.dropna())
    return lst

def plot_gc():
    #Variables that can be altered
    file_path = "GC_data/Woolly-mammoth_10000bpWS.txt"
    window_size = 10000 #This does not do anything except changes the plot title
    sample = 'Woolly mammoth' #This does not do anything except changes the plot title
    plot_color = '#1f7864'#'#f08b26'#'#0077b8
    numOfStds = 3

    data = readFiles(file_path)
    mean_GC = meanCalc(data)
    std_GC = stdCalc(data, mean_GC)
    findOutliers(data, mean_GC, std_GC, 1000000/window_size)

    removedOutliers = sum(data.loc[data['Outlier'], 'NumOfWindows'].tolist())
    totalWindows = data['NumOfWindows'].sum()
    removedOutliersFraction = 100*(removedOutliers/totalWindows)


    plt.subplot(2,1,1)
    plt.bar(data['GC%'], data['NumOfWindows'],color=plot_color, label='Included GC%')
    plt.bar(data['GC%'].where(data['Outlier']),data['NumOfWindows'],color='white', edgecolor='grey', label='Excluded GC%')
    plt.axvline(x=mean_GC, color='black')
    plt.axhline(y=1000000/window_size, color='grey')
    plt.text(mean_GC*1.01,max(data['NumOfWindows'])+max(data['NumOfWindows'])*0.001,f'Mean: {round(mean_GC,2)}')
    plt.text(max(data.index[data['Outlier']==False].tolist())*1.01,(1000000/window_size)*1.15,f'Min quantity threshold: {int(1000000/window_size)} \nWindows removed: {round(removedOutliersFraction,3)}% ({removedOutliers}/{totalWindows})')
    plt.xlabel('GC-content [%]')
    plt.ylabel('Quantity')
    plt.title(f'{sample}, window size {window_size} bp')
    plt.xlim([-5,105])
    plt.legend(fontsize=8)
    plt.yscale('log')


    plt.subplot(2,1,2)
    plt.errorbar(data['GC%'].where(data['Outlier'] == False),data['Mean'], yerr=data['Std'],color=plot_color, ecolor=plot_color, fmt='o',zorder=1)
    plt.errorbar(data['GC%'].where(data['Outlier']),data['Mean'], yerr=data['Std'], fmt='wo', ms=7, mec='grey',ecolor='grey',zorder=2)

    mean_list = dataFrameToList(data,'Mean')
    GC_list = dataFrameToList(data, 'GC%')
    # Fit trend line to data
    p = np.polyfit(GC_list,mean_list, deg=1) # Parameters from the fit of the polynomial
    y_model = np.polyval(p, GC_list) # Model the data using the parameters of the fitted straight line
    y_bar = np.mean(mean_list) # Mean
    R2 = np.sum((y_model - y_bar) ** 2) / np.sum((mean_list - y_bar) ** 2) # Coefficient of determination, R²
    xlim = plt.xlim()
    #Ploting the linear regression
    plt.plot(np.array(xlim), p[1] + p[0] * np.array(xlim), '--', linewidth=1.2, label=f'Line of Best Fit, R² = {R2:.2f}, y={round(p[0],4)}*x{round(p[1],4)}', color='black', zorder=3)
    plt.legend(fontsize=8)
    plt.axhline(y=0, color='black', linewidth=1)
    plt.xlim([-5,105])
    plt.xlabel('GC-content [%]')
    plt.ylabel('Comparative depth ratio')

    plt.show()


if __name__ == '__main__':
    plot_gc()