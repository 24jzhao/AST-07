#Idea.py: the work of a madman who wants to create 158 graphs at 4:38 AM
#Batch generate pdfs of light curves (~15min/curve)
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
#import astfunc -- TODO: add more stuff to it
#obj=4 -- unneeded here

rc('font', family='serif')
rc('mathtext', fontset='cm')
#full_file_list = os.listdir('../full') # creates a Python list containing the file paths for every object's time series

#A ******* disaster... but with finesse!


df = pd.read_csv('../csv/file_list.csv')
file_list = df.values.tolist()

df_l = pd.read_csv('../csv/lsa.csv') #Get objects to iterate over
lsa = df_l.values.tolist()

#make col 0 an int like it should be
for i in range(len(lsa)):
    del lsa[i][1]
    lsa[i] = int(lsa[i][0])

df_l1 = pd.read_csv('../csv/lsa1.csv') #Get master object indicies to find match names 
lsa1 = df_l1.values.tolist()
for i in range(len(lsa1)):
    lsa1[i][0] = int(lsa1[i][0])

#this only works when fnames is just the object name
#fnames = file_list[lsa1[lsa[obj]][0]]

df_cat = pd.read_csv('../csv/all_categ_vals.csv') #Get threes object incidies to match names 
all_categ_vals = df_cat.values.tolist()
all_categ_sort = [[w,x,y,z] for w,x,y,z in all_categ_vals if (x == 3 and y == 3 and z == 3)]

#Much cleaner way of using lsa, use when naming issue is solved
#fnames = (file_list[int(lsa[obj][0])])

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

def folded_light_curve(obj, period, pdf_name):
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
        plt.title(''+ " (" + str(obj) + ")") #TODO: make fnames[obj] work right UPDATE: it does work but I'm too lazy to fully implement it :P
        pdf.savefig()
        plt.close()

#outx = [] #Uncomment me if you are writing out a list of objects

periods = np.linspace(0.1, 1.0, 1000000) # This defines the search range of your period, you can specify it at your will. These are in days.
aval = 0 #aval = 61 last run

while aval < 158 #len(lsa), use int() if str or float 

#The bulk of it
    print(aval)
    obj = lsa[aval]
    t, mags, dy, filts = data_format(obj)
    LS_multi = periodic.LombScargleMultiband(Nterms_base=1, Nterms_band=0)#initiate structure variable
    LS_multi.fit(t, mags, dy, filts)#input our data
    P_multi = LS_multi.periodogram(periods)#function where input is periods
    plt.ion() #turns figure display on
    plt.figure()
    plt.scatter(periods, P_multi, s = 0.05)
    best_period = max(P_multi)
    for i in range(len(P_multi)):
        if P_multi[i] == best_period:
            index = i
    print(periods[index])
    #0.6 seconds ideal step
    # i/p_true + n where n is an integer all reciprocated is the beat frequency
    folded_light_curve(obj, periods[index], 'out/'+str(obj)+'.pdf')
    print('out/'+str(obj)+'.pdf generated.')
 

# Create file w/ indicies and object names
#    out1.append(aval)
#    outx.append([lsa[aval],lsa1[aval][0],file_list[lsa1[aval][0]]])
#    aval = aval+1
#pd.DataFrame(outx, columns=['3s_index','gen_index','name']).to_csv("outx.txt", index=False)
