
import pandas as pd
#from matplotlib import colors
import matplotlib.pyplot as plt
#import matplotlib.patches as patches
import numpy as np
from numpy.random import default_rng
from matplotlib.ticker import PercentFormatter

print('Done')

# |%%--%%| <XwDBIJQAig|SdkYwuYPQg>
# Basically a wrapper for a .csv containing a bunch of values

# df = pd.read_csv('../all_percentages_rand.csv') 
#NOTE: IN VIM $PWD is ~/, requires workaround:
df = pd.read_csv('ast07/AST-07/ast07/all_percentages_rand.csv')
all_percentages = df.values.tolist()
# Delete item 0 from every row because it messes stuff up
for i in range(len(all_percentages)):
    del(all_percentages[i][0])
    
print('Done')

# |%%--%%| <SdkYwuYPQg|hF2TfibX6j>

#Generate 10 random objects from the .csv
rng = default_rng() 
rint = rng.integers(low=0, high=26820, size=10)
# uncomment following line and add a 1d array of objects to reproduce a result:
#rint = ''#array goes here
# |%%--%%| <hF2TfibX6j|vOcVMzhnq9>

# Line 'em up...
for i in rint:
    print('Object '+str(i))
    fig, ax = plt.subplots()
    plt.grid(True)
    ax.set_axisbelow(True)
    ax.set_facecolor('#eaeaea') #BG Colo(u)r
    ax.set_autoscale_on(True) # Seems to do something good, it allegedly adjusts graph scope
    print(all_percentages[i])
    plt.xlabel('X-Label')
    plt.ylabel('Y-Label')

    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1)) # May or may not be accurate :/
    ax.set_title('RMS Percent for object '+str(i))
    ax.hist(all_percentages[i],bins=5,density='True',color='#4c72b0',rwidth=.7,histtype='barstacked',label='Inline label') # Addt'l options: ec='black',lw=2
    
    plt.show()

    

# |%%--%%| <vOcVMzhnq9|V2EjT7iSJx>

#... and knock em down!
#for i in range(len(all_percentages)):

for i in rint:
    labels = ['U', 'G', 'R', 'I', 'Z']
    x = 1+np.arange(len(labels)) #pos. of the blue guys
    fig, ax = plt.subplots()
    # Create bars
    barU = ax.bar(1, all_percentages[i][0], label='Test', color='purple')
    barG = plt.bar(2, all_percentages[i][1], label='Test', color='green')
    barR = plt.bar(3, all_percentages[i][2], label='Test', color='red')
    barI = plt.bar(4, all_percentages[i][3], label='Test', color='#4c72b0')
    barZ = plt.bar(5, all_percentages[i][4], label='Test', color='orange')
    #plt.bar(x, all_percentages[0], label='Test', color='#4c72b0')
    plt.grid(True)
    ax.set_axisbelow(True)
    ax.set_facecolor('#eaeaea') #BG Colo(u)r
    ax.set_autoscale_on(True) #Seems to do something good, it allegedly adjusts graph scope
    ax.set_xticks(x, labels)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1)) # May or may not be accurate :/
    ax.bar_label(barU, padding=1)
    ax.bar_label(barG, padding=1)
    ax.bar_label(barR, padding=1)
    ax.bar_label(barI, padding=1)
    ax.bar_label(barZ, padding=1)
    ax.set_title('RMS Percent for object '+str(i))

    

    plt.show()
            

# |%%--%%| <V2EjT7iSJx|U0Jxupz0WM>

df = pd.read_csv(‘all_categ_vals.csv’)
#df.drop(‘Unnamed: 0’,inplace=True,axis=1)
all_categ_vals = df.values.tolist()

# sorting algo!
all_categ_sort = [[w,x,y,z] for w,x,y,z in all_categ_vals if not (x==0 or y == 0 or z == 0)]

