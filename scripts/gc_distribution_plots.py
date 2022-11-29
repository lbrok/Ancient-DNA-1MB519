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

def stdCalc(data):
    sum_NumOfWindows = data.sum()['NumOfWindows']
    #Standard deviation calculations
    variance = 0
    for i in range(len(data)):
        variance += ((i-mean_GC)**2) * data['NumOfWindows'][i]
    variance = variance/sum_NumOfWindows
    std_GC = variance**0.5
    return std_GC

def findOutliers(data, mean, std):
    data['Outlier'] = ""
    for i in range(len(data)):
        if i < round(mean-numOfStds*std) or i > round(mean+numOfStds*std):
            data.loc[i,'Outlier'] = True
        else:
            data.loc[i,'Outlier'] = False

def dataFrameToList(data, column):
    lst = data[column].where(data['Outlier']==False)
    lst = list(lst.dropna())
    return lst

#Variables that can be altered
file_path = "GC_data\GC_GenomeWide_1000bpWS.txt"
#file_path = "GC_data\mammoth_1_1000bp.txt"
window_size = 1000 #This does not do anything except changes the plot title
sample = 'Pooled woolly mammoths' #This does not do anything except changes the plot title
plot_color = '#1f7864'#'#0077b8'
numOfStds = 3

data = readFiles(file_path)
mean_GC = meanCalc(data)
std_GC = stdCalc(data)
findOutliers(data, mean_GC, std_GC)

removedOutliers = sum(data.loc[data['Outlier'], 'NumOfWindows'].tolist())
totalWindows = data['NumOfWindows'].sum()
removedOutliersFraction = 100*(removedOutliers/totalWindows)

plt.subplot(2,1,1)
plt.bar(data['GC%'], data['NumOfWindows'],color=plot_color)
plt.bar(data['GC%'].where(data['Outlier']),data['NumOfWindows'],color='white', edgecolor='grey')
plt.axvline(x=mean_GC, color='black')
plt.axvline(x=round(mean_GC-numOfStds*std_GC)-0.5, color = 'grey')
plt.axvline(x=round(mean_GC+numOfStds*std_GC)+0.5, color = 'grey')
plt.text(mean_GC*1.01,max(data['NumOfWindows'])+max(data['NumOfWindows'])*0.001,f'Mean: {round(mean_GC,2)}')
plt.text((mean_GC+numOfStds*std_GC)*1.01,max(data['NumOfWindows'])*0.15,f'3 std confidence: {round(mean_GC,2)}±{round(numOfStds*std_GC,2)} GC% \nWindows removed: {round(removedOutliersFraction,2)}% ({removedOutliers}/{totalWindows})')
plt.xlabel('GC-content [%]')
plt.ylabel('Quantity')
plt.title(f'{sample}, window size {window_size} bp')
plt.xlim([-5,105])
#plt.yscale('log')

plt.subplot(2,1,2)
plt.errorbar(data['GC%'].where(data['Outlier'] == False),data['Mean'], yerr=data['Std'],color=plot_color, ecolor=plot_color, fmt='o',zorder=1)
plt.errorbar(data['GC%'].where(data['Outlier']),data['Mean'], yerr=data['Std'], fmt='wo', ms=7, mec='grey',ecolor='grey',zorder=2)
plt.axvline(x=round(mean_GC-numOfStds*std_GC)-0.5, color = 'grey')
plt.axvline(x=round(mean_GC+numOfStds*std_GC)+0.5, color = 'grey')
plt.text(mean_GC*1.01,max(data['NumOfWindows'])+max(data['NumOfWindows'])*0.005,f'Mean: {round(mean_GC,2)}')
plt.text((mean_GC+numOfStds*std_GC)*1.01,max(data['NumOfWindows'])*0.15,f'99% confidence: {round(mean_GC,2)}+-{round(numOfStds*std_GC,2)}')

mean_list = dataFrameToList(data,'Mean')
GC_list = dataFrameToList(data, 'GC%')
# Fit trend line to data
p = np.polyfit(GC_list,mean_list, deg=1) # Parameters from the fit of the polynomial
y_model = np.polyval(p, GC_list) # Model the data using the parameters of the fitted straight line
y_bar = np.mean(mean_list) # Mean
R2 = np.sum((y_model - y_bar) ** 2) / np.sum((mean_list - y_bar) ** 2) # Coefficient of determination, R²
xlim = plt.xlim()
#Ploting the linear regression
plt.plot(np.array(xlim), p[1] + p[0] * np.array(xlim), '--', linewidth=1.2, label=f'Line of Best Fit, R² = {R2:.2f}', color='black', zorder=3)
plt.legend(fontsize=8)
plt.axhline(y=0, color='black', linewidth=1)
plt.xlim([-5,105])
plt.xlabel('GC-content [%]')
plt.ylabel('Mammoth depth ratio - Elephant depth ratio / \nMammoth depth ratio + Elephant depth ratio')

plt.show()


