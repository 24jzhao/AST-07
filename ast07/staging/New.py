#Hopefully Generate working PDFs to out/
#Also test local libraries so clean code
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
from numpy.random import default_rng
#import locutils -- works
#import astfunc -- works
obj=4 #3s index no. 
rc('font', family='serif')
rc('mathtext', fontset='cm')
full_file_list = os.listdir('../full') # creates a Python list containing the file paths for every object's time series
#|%%--%%| <yckObVB2wW|CxZx0ewWbN>
#A ******* disaster...

df = pd.read_csv('../csv/file_list.csv')
file_list = df.values.tolist()
#retrieve list of known >20MOM vals
df_l = pd.read_csv('../csv/lsa.csv')
lsa = df_l.values.tolist()
#len(lsa)
#make col 0 an int like it should be
for i in range(len(lsa)):
    del lsa[i][1]
    lsa[i] = int(lsa[i][0])

df_l1 = pd.read_csv('../csv/lsa1.csv')
lsa1 = df_l1.values.tolist()
for i in range(len(lsa1)):
    lsa1[i][0] = int(lsa1[i][0])

#this only works when fnames is just the object name
#fnames = file_list[lsa1[lsa[obj]][0]]

print(lsa1)

df_cat = pd.read_csv('../csv/all_categ_vals.csv')
all_categ_vals = df_cat.values.tolist()
all_categ_sort = [[w,x,y,z] for w,x,y,z in all_categ_vals if (x == 3 and y == 3 and z == 3)]
#Much cleaner way of using lsa
#fnames = (file_list[int(lsa[obj][0])])

#|%%--%%| <CxZx0ewWbN|5eWx9TXGWn>
"""
mastlist=[]
for i in range(len(lsa)):
    for j in range(len(lsa1)):
        for k in range(len(file_list)):
            if lsa[i] == lsa1[i] 
"""
#|%%--%%| <5eWx9TXGWn|frvH1rjeDx>


all_mags_narrow = pd.read_csv('../csv/all_mags_narrow.csv')#.values.tolist()
all_mjd_narrow = pd.read_csv('../csv/all_mjd_narrow.csv')#.values.tolist()
all_magerrs_narrow = pd.read_csv('../csv/all_magerrs_narrow.csv')#.values.tolist()
all_mags_narrow.drop('Unnamed: 0',inplace=True,axis=1)
all_mjd_narrow.drop('Unnamed: 0',inplace=True,axis=1)
all_magerrs_narrow.drop('Unnamed: 0',inplace=True,axis=1)
all_mags = all_mags_narrow.values.tolist()
all_mjd = all_mjd_narrow.values.tolist()
all_magerrs = all_magerrs_narrow.values.tolist()
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
#|%%--%%| <frvH1rjeDx|ReShyMYaNw>
#generate random values to run light curves with
#rint = locutils.randomgen(lsa, 1)


#|%%--%%| <ReShyMYaNw|8yhsje2Zmb>
"""
def load_one_timeseries(file_name):
    timestamps = []
    mags = []
    magerrs = []
    expnums = []
    filters = []

    with open('../full/'+file_name) as data_file:
        first_line = data_file.readline()
        field = (first_line.split())[2]
        for line in data_file:
            if not line.startswith('#'):
                data = line.split()
                timestamps.append(data[0])
                mags.append(data[1])
                magerrs.append(data[2])
                expnums.append(data[3])
                filters.append(data[8])

    timestamps = np.asarray(timestamps, dtype=float)
    mags = np.asarray(mags, dtype=float)
    magerrs = np.asarray(magerrs, dtype=float)
    expnums = np.asarray(expnums, dtype=int)

    int_filters = np.empty(len(filters), dtype=int)
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=FutureWarning)
        int_filters[np.where(np.asarray(filters)=='U')[0]] = 0
        int_filters[np.where(np.asarray(filters)=='G')[0]] = 1
        int_filters[np.where(np.asarray(filters)=='R')[0]] = 2
        int_filters[np.where(np.asarray(filters)=='I')[0]] = 3
        int_filters[np.where(np.asarray(filters)=='I2')[0]] = 4
        int_filters[np.where(np.asarray(filters)=='Z')[0]] = 5

    num_entries = len(timestamps)

    return field, timestamps, mags, magerrs, expnums, int_filters, num_entries
"""

#|%%--%%| <8yhsje2Zmb|d439aYUHyM>
"""
def load_one_timeseries_field(file_name):
    timestamps = []
    mags = []
    magerrs = []
    expnums = []
    filters = []

    with open('../full/'+file_name) as data_file:
        first_line = data_file.readline()
        field = (first_line.split())[2]
        for line in data_file:
            if not line.startswith('#'):
                data = line.split()
                timestamps.append(data[0])
                mags.append(data[1])
                magerrs.append(data[2])
                expnums.append(data[3])
                filters.append(data[8])

    timestamps = np.asarray(timestamps, dtype=float)
    mags = np.asarray(mags, dtype=float)
    magerrs = np.asarray(magerrs, dtype=float)
    expnums = np.asarray(expnums, dtype=int)

    int_filters = np.empty(len(filters), dtype=int)
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=FutureWarning)
        int_filters[np.where(np.asarray(filters)=='U')[0]] = 0
        int_filters[np.where(np.asarray(filters)=='G')[0]] = 1
        int_filters[np.where(np.asarray(filters)=='R')[0]] = 2
        int_filters[np.where(np.asarray(filters)=='I')[0]] = 3
        int_filters[np.where(np.asarray(filters)=='I2')[0]] = 4
        int_filters[np.where(np.asarray(filters)=='Z')[0]] = 5

    num_entries = len(timestamps)

    return field, timestamps, mags, magerrs, expnums, int_filters, num_entries

"""
#|%%--%%| <d439aYUHyM|H4NDwEnVq2>
"""
def get_num_measurements(file_list, by_filter=False):
    if by_filter==False:
        num_measurements = np.empty(len(file_list))
        for i, fname in enumerate(file_list):
            _, _, _, _, _, _, num = load_one_timeseries(fname)
            num_measurements[i] = num
    elif by_filter==True:
        num_measurements = np.empty((len(file_list), 6))
        for i, fname in enumerate(file_list):
            _, _, _, _, _, filters, _ = load_one_timeseries(fname)
            for j in range(6):
                num_measurements[i,j] = np.sum(filters==j)
    return num_measurements


#|%%--%%| <H4NDwEnVq2|KISghI7R2g>
#df = pd.read_csv('../csv/file_list')
#file_list = df.tolist()

#if not file_list:
    #num_measurements = astfunc.get_num_measurements(full_file_list)
num_measurements = get_num_measurements(full_file_list)
file_list = (np.asarray(full_file_list))[num_measurements>15]
num_measurements_by_filter = get_num_measurements(file_list, by_filter=True)
    #num_measurements_by_filter = astfunc.get_num_measurements(file_list, by_filter=True)
file_list = file_list[np.all(num_measurements_by_filter>=15, axis=1)]
"""
#|%%--%%| <KISghI7R2g|1qr73Uiepz>

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

#|%%--%%| <1qr73Uiepz|VIRQ3CdD6W>
t, mags, dy, filts = data_format(obj)
obj = 0
print(obj)
periods = np.linspace(0.1, 1.0, 1000000) # This defines the search range of your period, you can specify it at your will. These are in days.

#|%%--%%| <VIRQ3CdD6W|dCi5Pk5Qef>

LS_multi = periodic.LombScargleMultiband(Nterms_base=1, Nterms_band=0)#initiate structure variable
LS_multi.fit(t, mags, dy, filts)#input our data
P_multi = LS_multi.periodogram(periods)#function where input is periods


#|%%--%%| <dCi5Pk5Qef|8Iu4cQ1yhl>

plt.ion() #turns figure display on
plt.figure()
plt.scatter(periods, P_multi, s = 0.05)
best_period = max(P_multi)
for i in range(len(P_multi)):
    if P_multi[i] == best_period:
        index = i
print(periods[index])


#|%%--%%| <8Iu4cQ1yhl|uI5MCQv74t>

#|%%--%%| <uI5MCQv74t|TLuSYeyvdV>


def folded_light_curve(obj, period, pdf_name):
    """
    Function to create folded light curve
    Parameters
    ---
    objects: int index of object of interest in file_list
    Returns
    ---
    pdf of magnitude vs MJD for all objects
    """
    plt.ioff()
    with PdfPages(pdf_name) as pdf:
        plt.figure(figsize = (9, 12))
        plt.xlabel('Phase')
        plt.ylabel('Magnitude')
        plt.gca().invert_yaxis()
        u_phase = [(i%period)/period for i in all_mjd[obj][0]]
        g_phase = [(i%period)/period for i in all_mjd[obj][1]]
        r_phase = [(i%period)/period for i in all_mjd[obj][2]]
        i_phase = [(i%period)/period for i in all_mjd[obj][3]]
        z_phase = [(i%period)/period for i in all_mjd[obj][4]]
        #f*p^2/length of observations, f is about 0.1, f = phase error for individual cycle
        plt.scatter(u_phase, all_mags[obj][0], s = 5, c = 'blue', label = 'u')
        plt.scatter(g_phase, all_mags[obj][1], s = 5, c = 'green', label = 'g')
        plt.scatter(r_phase, all_mags[obj][2], s = 5, c = 'purple', label = 'r')
        plt.scatter(i_phase, all_mags[obj][3], s = 5, c = 'gold', label = 'i')
        plt.scatter(z_phase, all_mags[obj][4], s = 5, c = 'tab:red', label = 'z')
        plt.legend()
        plt.title(''+ " (" + str(obj) + ")") #TODO: make fnames[obj] work right
        pdf.savefig()
        plt.close()
    #0.6 seconds ideal step
#|%%--%%| <TLuSYeyvdV|yv6qLpbnLI>
    #0.6 seconds ideal step
    # i/p_true + n where n is an integer all reciprocated is the beat frequency
#for i in lsa:
#        obj = lsa[i]
#        print(lsa[i])
folded_light_curve(obj, periods[index], 'out/'+str(obj)+'.pdf')
        print('out/'+str(obj)+'.pdf generated.')
