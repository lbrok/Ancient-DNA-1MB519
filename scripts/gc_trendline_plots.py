import matplotlib.pyplot as plt
from gc_distribution_plots import *

def create_vectors(k,m):
    X = [*range(-5,106)]
    Y = []
    for x in X:
        Y.append(k*x+m)
    return Y

def data_management(filePath, window_size=1000):
    data = readFiles(filePath)
    data_mean = meanCalc(data)
    data_std = stdCalc(data,data_mean)
    findOutliers(data, data_mean, data_std, 1000000/window_size)
    return data

def plot_scatter(x_vector, y_vector, data, color, trendline_label, transparency = 1):
    plt.plot(x_vector, y_vector, color, label=trendline_label, alpha=transparency)
    plt.scatter(data['GC%'].where(data['Outlier']),data['Mean'], c='white', alpha=transparency*0.5, edgecolors=color)
    plt.scatter(data['GC%'].where(data['Outlier'] == False), data['Mean'], c=color, alpha=transparency)


def plot_trendlines():
    neanderthal_color = '#f08b26'
    neanderthal_rev_color = '#3248a8'
    mammoth_color = '#1f7864'
    x_vector = [*range(-5,106)]
    neanderthal = create_vectors((-0.0035),0.1314)
    neanderthal_rev = create_vectors(0.0035,(-0.2151))
    mammoth = create_vectors(0.0066,(-0.2322))

    neanderthal_data = data_management("GC_data/Neanderthal_1000bpWS.txt", 1000)
    neanderthal_rev_data = data_management("GC_data/Neanderthal_1000bpWS_reversed.txt", 1000)
    mammoth_data = data_management("GC_data/Woolly-mammoth_1000bpWS.txt", 1000)

    plot_scatter(x_vector, neanderthal, neanderthal_data, neanderthal_color, 'Neanderthal, y=-0.0035x+0.1314', 0.5)
    plot_scatter(x_vector, mammoth, mammoth_data,mammoth_color,'Mammoth, y=0.0066x-0.2322')
    plot_scatter(x_vector, neanderthal_rev, neanderthal_rev_data, neanderthal_rev_color, 'Neanderthal reversed, y=0.0035x-0.2151')
    plt.xlim([-5,105])
    plt.ylim([-1,1])
    plt.axhline(y=0, color='black', linewidth = 1)
    plt.xlabel('GC-content [%]')
    plt.ylabel('Comparative depth ratio')
    plt.legend(fontsize=8)
    plt.show()

if __name__ == '__main__':
    plot_trendlines()