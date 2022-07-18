
import pandas as pd
import matplotlib.pyplot as plt
from pygnuplot import gnuplot
# |%%--%%| <XwDBIJQAig|SdkYwuYPQg>

df = pd.read_csv('~/ast07/all_percentages.csv')
all_percentages = df.values.tolist()

# |%%--%%| <SdkYwuYPQg|7DblBTbOJd>

def error_histogram(all_percentages): 
    for i in range(len(all_percentages)):
        plt.figure()
        plt.hist(all_percentages[i])
        plt.show()
        plt.savefig('test.png')
#|%%--%%| <7DblBTbOJd|vOcVMzhnq9>

