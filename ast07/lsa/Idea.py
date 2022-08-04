
#Idea.py: the work of a madman who wants to create 158 graphs at 4:38 AM
#Batch generate pdfs of light curves (~15min/curve)
#A ******* disaster... but with finesse!
#TODO: Add periodogram support

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
import a12lib

#dir w/ CSVs in relation to the notebook/.py (you need the slash at the end)
csvdir = '../csv/' 

#specify objects to look at here!
#obj0 = [18444, 25189, 21552, 12175, 24169, 149, 19875, 24483, 10753, 13369, 13729, 10753, 13369]

obj0 = [3933,15615,467,503,1041] #emergency revision of LSA PDFS
rc('font', family='sans')
rc('mathtext', fontset='cm')
df = pd.read_csv(csvdir+'file_list.csv')
file_list = df.values.tolist() #file list
#df_l = pd.read_csv(csvdir+'lsa.csv') #Get objects to iterate over
#lsa = df_l.values.tolist()
df_nd = pd.read_csv(csvdir+'all_threes_objects_filenames.csv')
df_nd = df_nd.drop(columns='Unnamed: 0')
valslist = df_nd.values.tolist()
oname = [] # object info (put it here to make it more readable)

####

#messy way to get fnames but it'll do
fnames = [] #get display name
#NOTE: to make this work using the full dataset, simply get a map of object indexs and names and replace valslist
for j in range(len(obj0)):
    for i in range(len(valslist)):
        if obj0[j] == valslist[i][1]:
            fnames.append(i)
## {{{ Onlyworks with three_object_indexes, but replacing 'em shouldn't hurt
all_mags_narrow = pd.read_csv(csvdir+'all_mags_narrow.csv')#.values.tolist()
all_mjd_narrow = pd.read_csv(csvdir+'all_mjd_narrow.csv')#.values.tolist()
all_magerrs_narrow = pd.read_csv(csvdir+'all_magerrs_narrow.csv')#.values.tolist()
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
## }}}

####

periods = np.linspace(0.1, 1.0, 1000000) # This defines the search range of your period, you can specify it at your will. These are in days.
print('[done] periods')

####

for a in range(len(obj0)):
    av1 = obj0[a] #object from user-input list
    an1 = fnames[a] #index of the object
    obj = an1 #make the code work
    for i in range(len(valslist)):
        if av1 == valslist[i][1]:
            oname = (i,valslist[i][0]) #temp list containing relative index and name (for display)
#    break #uncomment this if testing something, DO NOT split cells
#The bulk of it
    t, mags, dy, filts = a12lib.data_format(obj, all_mags, all_mjd, all_magerrs)
    print('[done] data_format')
    LS_multi = periodic.LombScargleMultiband(Nterms_base=1, Nterms_band=0)#initiate structure variable
    LS_multi.fit(t, mags, dy, filts)#input our data
    P_multi = LS_multi.periodogram(periods)#function where input is periods
    print('[done] multi-stuff')
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
    a12lib.folded_light_curve(obj, periods[index], 'out/num/gen/'+str(av1)+'.pdf', oname, av1, all_mjd, all_mags)
    a12lib.folded_light_curve(obj, periods[index], 'out/nam/gen/'+str(oname[1])+'.pdf', oname, av1, all_mjd, all_mags)
    print(str(oname[1])+'.pdf generated.')
    print(str(av1)+'.pdf generated.')
print(oname)
