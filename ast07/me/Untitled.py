
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

df = pd.read_csv('../all_percentages_rand.csv') 
#NOTE: IN VIM $PWD is ~/, requires workaround:
#df = pd.read_csv('ast07/AST-07/ast07/all_percentages_rand.csv')
all_percentages = df.values.tolist()

# Delete item 0 from every row because it messes stuff up
for i in range(len(all_percentages)):
    del(all_percentages[i][0])
    del(all_percentages[i][3]) # IZ unused, so deleting them (works fine)
    del(all_percentages[i][3])
    
#    all_percentages[i][3] = 0 # Alternate way to delete them without breaking anything
#    all_percentages[i][4] = 0

print(all_percentages[2740])
print('Done')

# |%%--%%| <SdkYwuYPQg|hF2TfibX6j>

#Generate 10 random objects from the .csv
rng = default_rng() 
rint = rng.integers(low=0, high=26820, size=10)
print('Objects Observed: '+str(rint))
# uncomment following line and add a 1d array of objects to reproduce a result:
#rint = ''#array goes here

# |%%--%%| <hF2TfibX6j|vOcVMzhnq9>

# Unused because it is inaccurate, left in for compatability

# Quick explanation why this looks inaccurate:
# Histograms graph a range of values (bin width) against frequency;
# So the Y-Vals in the graph is the number of occurances the range of X-values appear in the dataset (printed above graph)
# Odd thing that should be looked into: sum of Y-Vals always equals 3

# Line em up...
#for i in rint:
#    print('Object '+str(i))
#    fig, ax = plt.subplots()
#    plt.grid(True)
#    ax.set_axisbelow(True)
#    ax.set_facecolor('#eaeaea') #BG Colo(u)r
#    ax.set_autoscale_on(True) # Seems to do something good, it allegedly adjusts graph scope
#    print(all_percentages[i])
#    plt.xlabel('X-Label')
#    plt.ylabel('Y-Label')
#
#    ax.set_title('RMS Error Percent for object '+str(i))
#    ax.hist(all_percentages[i],bins=len(all_percentages[i]),color='#4c72b0',rwidth=.7,histtype='barstacked',label='Inline label') # Addt'l options: ec='black',lw=2
    
#    plt.show()
print('WARN: This feature is deprecated and may be removed in the future')
    

# |%%--%%| <vOcVMzhnq9|V2EjT7iSJx>

#... and knock em down!
# Error2Bar

#for i in range(len(all_percentages)):

for i in rint:
#    labels = ['U', 'G', 'R', 'I', 'Z'] outdated label list
    labels = ['U', 'G', 'R']
    x = 1+np.arange(len(labels)) #pos. of the blue guys
    fig, ax = plt.subplots()
    # Create bars (ax.bar and plt.bar are interchangable as long as ax is defined)
    barU = ax.bar(1, all_percentages[i][0], label='Test', color='purple')
    barG = plt.bar(2, all_percentages[i][1], label='Test', color='green')
    barR = plt.bar(3, all_percentages[i][2], label='Test', color='red')
#    barI = plt.bar(4, all_percentages[i][3], label='Test', color='#4c72b0') # Kept incase IZ needs to be revisited
#    barZ = plt.bar(5, all_percentages[i][4], label='Test', color='orange')

    plt.grid(True)
    ax.set_axisbelow(True)
    ax.set_facecolor('#eaeaea') #BG Colo(u)r
    ax.set_autoscale_on(True) #Seems to do something good, it allegedly adjusts graph scope
    ax.set_xticks(x, labels)
  #  ax.yaxis.set_major_formatter() # Figure out how to make it say '%' at the end
    ax.bar_label(barU, padding=1)
    ax.bar_label(barG, padding=1)
    ax.bar_label(barR, padding=1)
#    ax.bar_label(barI, padding=1) # See barI and barZ notes
#    ax.bar_label(barZ, padding=1) # 
    ax.set_title('RMS Percent Error for object '+str(i))

    

    plt.show()
            

# |%%--%%| <V2EjT7iSJx|U0Jxupz0WM>

df = pd.read_csv('../all_categ_vals.csv')
#df.drop('Unnamed: 0',inplace=True,axis=1) # Don't cut off object numbers
all_categ_vals = df.values.tolist()

# sorting algo!
# NOTE: keeps object number as item 0, chop off 'w' if you don't want object no.
all_categ_sort = [[w,x,y,z] for w,x,y,z in all_categ_vals if (x==3 and y==3 and z ==3)]


# |%%--%%| <U0Jxupz0WM|BnCi5as2ao>

#for i in range(len(all_categ_sort)):
#    print(all_categ_sort[i][0])

# |%%--%%| <BnCi5as2ao|pSqLIzg7zg>


