"""°°°
# OVERVIEW:
The data required for this project can be found at this link: 
https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/cfhtls/dfspt.html
One must download the DeepVarFull.tar.gz file and drag it/download it into the Jupyter notebook or project to work with it. The data should contain approximately 28000 files to be extracted in which each file represents an astronomical object. There is information about the object within each file, like its magnitude measurements and the filters in which the measurements were taken, etc.. Each of the astronomical objects had several brightness measurements taken in each of the six filters ('U','G','R','I1','I2', and 'Z' are our names for them).
In the code below, you will find that we have combined the fourth and fifth filters (i1 and i2) into a single filter, thus having more magnitude measurements than the other filters.
°°°"""
# |%%--%%| <8hwXvaNZnh|NjtQqynrtp>

#downloading gatsby

# |%%--%%| <NjtQqynrtp|8X5ZoLSZ5N>

import sys
#import statements
import numpy as np
import math
import pandas as pd
import os
import warnings
from astropy.stats import sigma_clip
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy.ma as ma
import copy
from matplotlib import rc
import numpy.ma as ma
import matplotlib as mpl
import sklearn.metrics
from itertools import chain
from gatspy import datasets, periodic
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LinearRegression

rc('font', family='serif')
rc('mathtext', fontset='cm')

full_file_list = os.listdir('../full') # creates a Python list containing the file paths for every object's time series

# |%%--%%| <8X5ZoLSZ5N|YN09D7F4yY>
# Cool way to sort array in only one line
df = pd.read_csv('csv/all_categ_vals.csv')
all_categ_vals = df.values.tolist()

#Get indicies for value set
three_object_indexes = [w for w,x,y,z in all_categ_vals if (x==3 and y == 3 and z == 3)]

# |%%--%%| <YN09D7F4yY|I219QxkTi5>

all_mags_narrow = pd.read_csv('csv/all_mags_narrow.csv')#.values.tolist()
all_mjd_narrow = pd.read_csv('csv/all_mjd_narrow.csv')#.values.tolist()
all_magerrs_narrow = pd.read_csv('csv/all_magerrs_narrow.csv')#.values.tolist()
all_mags_narrow.drop('Unnamed: 0',inplace=True,axis=1)
all_mjd_narrow.drop('Unnamed: 0',inplace=True,axis=1)
all_magerrs_narrow.drop('Unnamed: 0',inplace=True,axis=1)
all_mags = all_mags_narrow.values.tolist()
all_mjd = all_mjd_narrow.values.tolist()
all_magerrs = all_magerrs_narrow.values.tolist()
#|%%--%%| <I219QxkTi5|HGWxOmHzRJ>

for i in range(len(all_mags)):
    for j in range(len(all_mags[i])):
        all_mags[i][j] = all_mags[i][j][1:len(all_mags[i][j])-1].replace('\n','').split()
        all_mjd[i][j] = all_mjd[i][j][1:len(all_mjd[i][j])-1].replace('\n','').split()
        all_magerrs[i][j] = all_magerrs[i][j][1:len(all_magerrs[i][j])-1].replace('\n','').split()
for i in range(len(all_mags)):
    for j in range(len(all_mags[i])):
        for z in range(len(all_mags[i][j])):
            all_mags[i][j][z] = float(all_mags[i][j][z])
            all_mjd[i][j][z] = float(all_mjd[i][j][z])
            all_magerrs[i][j][z] = float(all_magerrs[i][j][z])


#|%%--%%| <HGWxOmHzRJ|exHE8NGZMb>

def data_format (object_number):
    filts = []
    temp_u = ['u' for i in all_mags[object_number][0]]
    temp_g = ['g' for i in all_mags[object_number][1]]
    temp_r = ['r' for i in all_mags[object_number][2]]
    temp_i = ['i' for i in all_mags[object_number][3]]
    temp_z = ['z' for i in all_mags[object_number][4]]
    filts.append(temp_u)
    filts.append(temp_g)
    filts.append(temp_r)
    filts.append(temp_i)
    filts.append(temp_z)
    filts = list(chain.from_iterable(filts))
    mags = list(chain.from_iterable(all_mags[object_number])) #turning 2d array into 1d
    t = list(chain.from_iterable(all_mjd[object_number]))
    #dy = [1 for i in mags]
    dy = list(chain.from_iterable(all_magerrs[object_number]))
    dy = [99 if i == 0 else i for i in dy]
    return t, mags, dy, filts

#|%%--%%| <exHE8NGZMb|XUbQvcqmNO>

def lomb_scargle_list (object_list):
    print(i) #not sure if this actually works, but oh well
    max_over_mean = []
    for cur_object in object_list:
        t, mags, dy, filts = data_format(cur_object)
        periods = np.linspace(0.1, 1.0, 100000)
        LS_multi = periodic.LombScargleMultiband(Nterms_base=1, Nterms_band=0)#initiate structure variable
        LS_multi.fit(t, mags, dy, filts)#input our data
        P_multi = LS_multi.periodogram(periods)#function where input is periods
        max_over_mean.append(np.max(P_multi)/np.mean(P_multi))
    return max_over_mean,object_list

#|%%--%%| <XUbQvcqmNO|xSkcEvJdom>
print('lsa')
max_over_mean_threes,object_list = lomb_scargle_list([i for i in range(397)])
print('lsa done')
print(object_list_saved)
max_over_mean_threes_saved = max_over_mean_threes
object_list_saved = object_list
#|%%--%%| <xSkcEvJdom|ZKXubZNHmX>
#print(max_over_mean_threes)
narrowed_objects = []
for i in range(len(max_over_mean_threes_saved)):
    if max_over_mean_threes_saved[i] >= 20:
        narrowed_objects.append([object_list_saved[i],max_over_mean_threes_saved[i]])


len(narrowed_objects)

#|%%--%%| <ZKXubZNHmX|5hvgRVGBTX>
import csv
#dfx = pd.DataFrame(narrowed_objects)
with open('blahblah', 'w', newline='') as blah:
    wr = csv.writer(blah, quoting=csv.QUOTE_ALL)
    wr.writerow(narrowed_objects)


#|%%--%%| <5hvgRVGBTX|GcNKluNsHb>
df_x = pd.DataFrame(narrowed_objects)
import base64
import pandas as pd
from IPython.display import HTML
def create_download_link( df, title = "Download CSV file", filename = "data.csv"):
    csv = df.to_csv()
    b64 = base64.b64encode(csv.encode())
    payload = b64.decode()
    html = '<a download="{filename}" href="data:text/csv;base64,{payload}" target="_blank">{title}</a>'
    html = html.format(payload=payload,title=title,filename=filename)
    return HTML(html)
create_download_link(df_x)

