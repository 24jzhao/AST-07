#!/usr/bin/env python
# coding: utf-8

# # OVERVIEW:
# The data required for this project can be found at this link: 
# https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/cfhtls/dfspt.html
# One must download the DeepVarFull.tar.gz file and drag it/download it into the Jupyter notebook or project to work with it. The data should contain approximately 28000 files to be extracted in which each file represents an astronomical object. There is information about the object within each file, like its magnitude measurements and the filters in which the measurements were taken, etc.. Each of the astronomical objects had several brightness measurements taken in each of the six filters ('U','G','R','I1','I2', and 'Z' are our names for them).
# In the code below, you will find that we have combined the fourth and fifth filters (i1 and i2) into a single filter, thus having more magnitude measurements than the other filters.

# In[1]:


#downloading gatsby
import sys
get_ipython().system('{sys.executable} -m pip install gatspy')


# In[2]:


#fields[12]


# In[3]:


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









#get_ipython().run_line_magic('matplotlib', 'notebook')

full_file_list = os.listdir('../full') # creates a Python list containing the file paths for every object's time series




# The following code is the basis for all the analysis.

# In[4]:


def load_one_timeseries(file_name):
    """
    Function to access light curve data for one object.
    
    Parameters
    ---
    file_name : str
        File path to CFHTLS time series data file (eg. 'CFHTLS-VAR-J022359.77-041810.1.mjdmag')
        
    Returns
    ---
    field : str
        CFHTLS field (D1, D2, D3, or D4) object is found in
    timestamps : ndarray
        Array of recorded times for each measurement
    mags : ndarray
        Array of measured magnitude for each measurement
    magerrs : ndarray
        Array of measurement error in magnitude for each measurement
    expnums : ndarray
        Array of the MegaCam exposure number in which each measurement can be found
    int_filters : ndarray
        Array of integers representing each CFHTLS filter indicating what filter each magnitude was measured in
        (u = 0, g = 1, r = 2, i1 = 3, i2 = 4, z = 5 - there were two different I-band filters)
    num_entries : int
        Total number of measurements
    
    """
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


# In[5]:


def load_one_timeseries_field(file_name):
    """
    Function to access light curve data for one object.
    
    Parameters
    ---
    file_name : str
        File path to CFHTLS time series data file (eg. 'CFHTLS-VAR-J022359.77-041810.1.mjdmag')
        
    Returns
    ---
    field : str
        CFHTLS field (D1, D2, D3, or D4) object is found in
    timestamps : ndarray
        Array of recorded times for each measurement
    mags : ndarray
        Array of measured magnitude for each measurement
    magerrs : ndarray
        Array of measurement error in magnitude for each measurement
    expnums : ndarray
        Array of the MegaCam exposure number in which each measurement can be found
    int_filters : ndarray
        Array of integers representing each CFHTLS filter indicating what filter each magnitude was measured in
        (u = 0, g = 1, r = 2, i1 = 3, i2 = 4, z = 5 - there were two different I-band filters)
    num_entries : int
        Total number of measurements
    
    """
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


# In[6]:


def get_num_measurements(file_list, by_filter=False):
    """
    Function to obtain the number of measurements in each time series file.
    
    Parameters
    ---
    file_list : list of str
        List of file paths (eg. one created by os.listdir(folder_path))
    by_filter : bool
        Defaults to False. Set to true to return an array containing number of measurements per filter
        
    Returns
    ---
    num_measurements : ndarray
        Array containing number of time series measurements per object
        If by_filter=False, Nx1 array, where N = number of objects
        If by_filter=True, NxM array, where N = number of objects and M = number of filters
    
    """
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


# In[7]:


num_measurements = get_num_measurements(full_file_list)


# In[8]:


file_list = (np.asarray(full_file_list))[num_measurements>15]


# In[9]:


num_measurements_by_filter = get_num_measurements(file_list, by_filter=True)


# In[10]:


#there must be at least 15 measurements per band for an object for us to use it
file_list = file_list[np.all(num_measurements_by_filter>=15, axis=1)]


# In[11]:




def get_rms_by_median(magnitudes, bool_array):
    true_mags =  [mag_val for (mag_val,bool_val) in zip(magnitdues,bool_array) if bool_val == True]
    median = np.median(true_magnitudes)
    return math.sqrt(np.mean([((i - median) * (i-median) for i in true_mags)]))
    
    
def convert_i1_to_i2(i1_magnitudes): #converts i1 magnitudes to i2
                     
    i2_magnitudes = []
    x = range(len(i1_magnitudes))
    for i in x:
        conversion = 0.04 * (medians[i, 1] - medians[i, 3]) + i1_magnitudes[i] - 0.076
        i2_magnitudes.append(conversion)
    return i2_magnitudes  
    
  
def get_medians(file_list):
    medians = np.empty((len(file_list), 6))
    for i, fname in enumerate(file_list):
        _, _, mags, _, _, filters, _ = load_one_timeseries(fname)
        for j in range(6):
            medians[i,j] = np.median(mags[filters==j])
    return medians

def get_boolean_array_5 (boolean_array, magnitudes,sigma_val):
    """
    Function that returns boolean values indicating whether measurement was sigma clipped (true = not sigma clipped; false = sigma clipped)
    Parameters
    ---
    boolean_array: array of true false values (i.e.: u_boolean, r_boolean)
    magnitudes: array of magnitudes (i.e.: u_mags, r_mags)
    Returns
    ---
    new_boolean_array: array of true false values
    """
    true_mags = []
    true_mags = [boolean_array[i] for i in range(len(magnitudes)) if boolean_array[i] == True]
    rms = get_rms_by_median(magnitudes, boolean_array)
    for i in range (len(magnitudes)):
        if boolean_array[i] == True:
            if abs((magnitudes[i] - np.median(true_mags))/rms) > sigma_val:
                boolean_array[i] = False
    return boolean_array

                     
medians = get_medians(file_list)

    


# In[12]:


#downloading data
magnitudes_raw_data = []
filters_raw_data = []
fnames = []
fields = []
mjds = []
magerrs = []
for fname in file_list:
    field, mjd, magnitudes, magerr, _, filters, _ = load_one_timeseries(fname)
    magnitudes_raw_data.append(magnitudes)
    filters_raw_data.append(filters)
    fnames.append(fname)
    fields.append(field)
    mjds.append(mjd)
    magerrs.append(magerr)
    
    


# In[13]:


#creating 6 empty lists for each object, each object is an element in the overall 3d list
all_mags = []
for i in range(len(magnitudes_raw_data)):
    all_mags.append([[],[],[],[],[],[]])
    
all_mjd = []
for i in range(len(mjds)):
    all_mjd.append([[],[],[],[],[],[]])
    
all_magerrs = []
for i in range(len(mjds)):
    all_magerrs.append([[],[],[],[],[],[]])



# In[14]:


#deleting '.mjdmag' from file name of each object
for i, fname in enumerate(fnames):
    fnames[i] = fname.replace(".mjdmag", "")


# In[15]:


for i in range(len(fnames)):
    if fnames[i] == "CFHTLS-VAR-J221527.94-180359.8":
        print(i)


# In[16]:


#sort measurements by filter
for i in range(len(magnitudes_raw_data)):
    all_mags[i][0] = [magnitudes_raw_data[i][j] for j in range(len(magnitudes_raw_data[i])) if filters_raw_data[i][j] == 0] #u
    all_mags[i][1] = [magnitudes_raw_data[i][j] for j in range(len(magnitudes_raw_data[i])) if filters_raw_data[i][j] == 1] #g
    all_mags[i][2] = [magnitudes_raw_data[i][j] for j in range(len(magnitudes_raw_data[i])) if filters_raw_data[i][j] == 2] #r
    all_mags[i][3] = [magnitudes_raw_data[i][j] for j in range(len(magnitudes_raw_data[i])) if filters_raw_data[i][j] == 3] #i1
    all_mags[i][4] = [magnitudes_raw_data[i][j] for j in range(len(magnitudes_raw_data[i])) if filters_raw_data[i][j] == 4] #i2
    all_mags[i][5] = [magnitudes_raw_data[i][j] for j in range(len(magnitudes_raw_data[i])) if filters_raw_data[i][j] == 5] #z


#i_one_converted = convert_i1_to_i2(i_one_mags)
#i_two_mags.extend(i_one_converted)

for i in range(len(mjds)):
    all_mjd[i][0] = [mjds[i][j] for j in range(len(mjds[i])) if filters_raw_data[i][j] == 0]
    all_mjd[i][1] = [mjds[i][j] for j in range(len(mjds[i])) if filters_raw_data[i][j] == 1]
    all_mjd[i][2] = [mjds[i][j] for j in range(len(mjds[i])) if filters_raw_data[i][j] == 2]
    all_mjd[i][3] = [mjds[i][j] for j in range(len(mjds[i])) if filters_raw_data[i][j] == 3]
    all_mjd[i][4] = [mjds[i][j] for j in range(len(mjds[i])) if filters_raw_data[i][j] == 4]
    all_mjd[i][5] = [mjds[i][j] for j in range(len(mjds[i])) if filters_raw_data[i][j] == 5]
    
for i in range(len(magerrs)):
    all_magerrs[i][0] = [magerrs[i][j] for j in range(len(magerrs[i])) if filters_raw_data[i][j] == 0]
    all_magerrs[i][1] = [magerrs[i][j] for j in range(len(magerrs[i])) if filters_raw_data[i][j] == 1]
    all_magerrs[i][2] = [magerrs[i][j] for j in range(len(magerrs[i])) if filters_raw_data[i][j] == 2]
    all_magerrs[i][3] = [magerrs[i][j] for j in range(len(magerrs[i])) if filters_raw_data[i][j] == 3]
    all_magerrs[i][4] = [magerrs[i][j] for j in range(len(magerrs[i])) if filters_raw_data[i][j] == 4]
    all_magerrs[i][5] = [magerrs[i][j] for j in range(len(magerrs[i])) if filters_raw_data[i][j] == 5]


# In[17]:


#convert i1 to i2 through predefined function 'convert_i1_to_i2'
for i in range(len(all_mags)):
    i_one_converted = convert_i1_to_i2(all_mags[i][3])
    all_mags[i][4].extend(i_one_converted)
    
for i in range(len(all_mjd)):
    i_one_mjd = all_mjd[i][3]
    all_mjd[i][4].extend(i_one_mjd)
    
for i in range(len(all_magerrs)):
    i_one_magerr = all_magerrs[i][3]
    all_magerrs[i][4].extend(i_one_magerr)


# In[18]:


#remove original i1 measurements from 3d list
for i in range(len(all_mags)):
    all_mags[i].remove(all_mags[i][3])
    
for i in range(len(all_mjd)):
    all_mjd[i].remove(all_mjd[i][3])
    
for i in range(len(all_magerrs)):
    all_magerrs[i].remove(all_magerrs[i][3])


# In[19]:


#create a copy of all_mags
all_mags_2 = []
for i in range(len(all_mags)):
    all_mags_2.append([[],[],[],[],[]])


for i in range(len(all_mags)):
    for j in range(len(all_mags[i])):
            all_mags_2[i][j] = copy.deepcopy(all_mags[i][j])


# In[20]:


all_mags_copy = copy.deepcopy(all_mags)


# In[21]:


#5 sigma clip all_mags 
for i in range(len(all_mags)):
    all_mags[i][0] = sigma_clip(all_mags[i][0],sigma=4,maxiters=3,masked=False,copy=False)
    all_mags[i][1] = sigma_clip(all_mags[i][1],sigma=4,maxiters=3,masked=False,copy=False)
    all_mags[i][2] = sigma_clip(all_mags[i][2],sigma=4,maxiters=3,masked=False,copy=False)
    all_mags[i][3] = sigma_clip(all_mags[i][3],sigma=4,maxiters=3,masked=False,copy=False)
    all_mags[i][4] = sigma_clip(all_mags[i][4],sigma=4,maxiters=3,masked=False,copy=False)
    


# In[22]:


#5 sigma clip all_mags_copy and RETURN MASKED
for i in range(len(all_mags_copy)):
    all_mags_copy[i][0] = sigma_clip(all_mags_copy[i][0],sigma=4,maxiters=3,masked=True,copy=False)
    all_mags_copy[i][1] = sigma_clip(all_mags_copy[i][1],sigma=4,maxiters=3,masked=True,copy=False)
    all_mags_copy[i][2] = sigma_clip(all_mags_copy[i][2],sigma=4,maxiters=3,masked=True,copy=False)
    all_mags_copy[i][3] = sigma_clip(all_mags_copy[i][3],sigma=4,maxiters=3,masked=True,copy=False)
    all_mags_copy[i][4] = sigma_clip(all_mags_copy[i][4],sigma=4,maxiters=3,masked=True,copy=False)
    


# In[23]:


#deleting corresponding sigma-clipped out data in all_mjd and all_magerrs
for i in range(len(all_mags_copy)):
    for j in range(len(all_mags_copy[i])):
        temp_mask = ma.getmaskarray(all_mags_copy[i][j])
        all_mjd[i][j] = np.delete(all_mjd[i][j], temp_mask)
        all_magerrs[i][j] = np.delete(all_magerrs[i][j], temp_mask)


# In[24]:


#4.5 sigma clip all_mags_2
for i in range(len(all_mags_2)):
    all_mags_2[i][0] = sigma_clip(all_mags_2[i][0],sigma=3,maxiters=3,masked=False,copy=False)
    all_mags_2[i][1] = sigma_clip(all_mags_2[i][1],sigma=3,maxiters=3,masked=False,copy=False)
    all_mags_2[i][2] = sigma_clip(all_mags_2[i][2],sigma=3,maxiters=3,masked=False,copy=False)
    all_mags_2[i][3] = sigma_clip(all_mags_2[i][3],sigma=3,maxiters=3,masked=False,copy=False)
    all_mags_2[i][4] = sigma_clip(all_mags_2[i][4],sigma=3,maxiters=3,masked=False,copy=False)


# In[25]:


#separating objects by field
#these lists have indicies of objects
fields_1 = [i for i in range(len(fields)) if fields[i] == "D1"]
fields_2 = [i for i in range(len(fields)) if fields[i] == "D2"]
fields_3 = [i for i in range(len(fields)) if fields[i] == "D3"]
fields_4 = [i for i in range(len(fields)) if fields[i] == "D4"]


# In[26]:


#separating measurements of objects by fields
D1 = [all_mags[i] for i in range(len(all_mags)) if i in fields_1]
D2 = [all_mags[i] for i in range(len(all_mags)) if i in fields_2]
D3 = [all_mags[i] for i in range(len(all_mags)) if i in fields_3]
D4 = [all_mags[i] for i in range(len(all_mags)) if i in fields_4]


# In[27]:


#separating 4.5 sigma clipped measurements of objects by fields
D1_sigma = [all_mags_2[i] for i in range(len(all_mags_2)) if i in fields_1]
D2_sigma = [all_mags_2[i] for i in range(len(all_mags_2)) if i in fields_2]
D3_sigma = [all_mags_2[i] for i in range(len(all_mags_2)) if i in fields_3]
D4_sigma = [all_mags_2[i] for i in range(len(all_mags_2)) if i in fields_4]


# In[28]:


#sorting file names by field
D1_fnames = [fnames[i] for i in range(len(fnames)) if i in fields_1]
D2_fnames = [fnames[i] for i in range(len(fnames)) if i in fields_2]
D3_fnames = [fnames[i] for i in range(len(fnames)) if i in fields_3]
D4_fnames = [fnames[i] for i in range(len(fnames)) if i in fields_4]


# In[47]:


def get_categ_vals(band, mags, mags_for_median):
    """
    Function that plots median magnitude vs. decimal log rms for a given filter and given field.
    Parameters
    ---
    band: int for filter (0 = u, 1 = g, ..., 4 = z)
    mags: choose from D1, D2, D3, D4
    mags_for_median: choose from 4.5 sigma clipped mags; choose from D1_sigma, D2_sigma, D3_sigma, D4_sigma
    Returns
    ---
    Plot of median magnitudes vs. log rms for a specific band and specific field
    """
    all_true_mags = []
    median_mags = []
    medians = []
    all_rms = []
    all_rms_two = []
    medians_for_pink = []
    """
    for i in range(len(mags)):
        true_mags = [mags[i][band][j] for j in range(len(mags[i][band]))]
        all_true_mags.append(true_mags)
    for i in range(len(mags_for_median)):
        temp_mags_for_median = [mags_for_median[i][band][j] for j in range(len(mags_for_median[i][band]))]
        median_mags.append(temp_mags_for_median)
    for i in range(len(all_true_mags)):
        squared_differences = []
        median = np.median(all_true_mags[i])
        medians.append(median)
        for j in range(len(all_true_mags[i])):
            difference = all_true_mags[i][j] - median
            squared_differences.append(difference * difference)
            rms = np.mean(squared_differences) ** 0.5
        all_rms.append(rms)
    for i in range(len(median_mags)):
        squared_differences_two = []
        median_two = np.median(median_mags[i])
        medians_for_pink.append(median_two)
        for j in range(len(median_mags[i])):
            difference = median_mags[i][j] - median_two
            squared_differences_two.append(difference * difference)
            rms = np.mean(squared_differences_two) ** 0.5
        all_rms_two.append(rms)
    """
    #simplified
    for i in range(len(mags)):
        true_mags = [mags[i][band][j] for j in range(len(mags[i][band]))]
        all_true_mags.append(true_mags)
        temp_mags_for_median = [mags_for_median[i][band][j] for j in range(len(mags_for_median[i][band]))]
        median_mags.append(temp_mags_for_median)
    ####
    for i in range(len(all_true_mags)): #going through 5-sigma clipped mags (i is each object)
        squared_differences = []
        median = np.median(all_true_mags[i]) #median for 1 object's magnitudes in a specific filter (see: band parameter)
        medians.append(median) #medians = x-vals for blue dots
        differences = [all_true_mags[i][j] - median for j in range(len(all_true_mags[i]))]
        squared_differences.extend([i * i for i in differences])
        rms = np.mean(squared_differences) ** 0.5
        all_rms.append(rms)
        #second part
        squared_differences_two = []
        median_two = np.median(median_mags[i])
        medians_for_pink.append(median_two)
        differences_two = [median_mags[i][j] - median_two for j in range(len(median_mags[i]))]
        squared_differences_two.extend([i * i for i in differences_two])
        rms = np.mean(squared_differences_two) ** 0.5
        all_rms_two.append(rms)
    
    red_rms = []
    in_bin_medians = []
    bounds = np.arange(min(medians_for_pink), max(medians_for_pink), 0.2)
    for i in range(len(bounds) - 1):
        in_bin = []
        temp_rms = []
        for j in range(len(medians_for_pink)):
            if(bounds[i] < medians_for_pink[j] and bounds[i+1] > medians_for_pink[j]):
                in_bin.append(medians_for_pink[j])
                temp_rms.append(all_rms_two[j])
        in_bin_medians.append(np.median(in_bin))
        red_rms.append(np.median(temp_rms))
    nan_bool_x = np.isnan(in_bin_medians)
    nan_bool_y = np.isnan(red_rms)
    in_bin_medians = np.delete(in_bin_medians, nan_bool_x)
    red_rms = np.delete(red_rms, nan_bool_y)
    for_color = np.empty_like(red_rms)
    plt.figure(figsize = (9, 12))
    plt.scatter(medians,np.log10(all_rms), s = 0.20)
    plt.scatter(in_bin_medians,np.log10(red_rms), s = 13, c = "crimson")#c = for_color, cmap = 'PuRd'
    if band == 0:
        plt.ylabel(r'$\sigma_%s$' % 'u')
        plt.xlabel(r'$<%s>$' % 'u')
        if mags == D1:
            X = (in_bin_medians[6:len(in_bin_medians)-3])
            y = np.log10(red_rms[6:len(red_rms) - 3])
        elif mags == D2:
            X = (in_bin_medians[6:len(in_bin_medians)-5])
            y = np.log10(red_rms[6:len(red_rms) - 5])
        elif mags == D3:
            X = (in_bin_medians[10:len(in_bin_medians)-6])
            y = np.log10(red_rms[10:len(red_rms) - 6])
        elif mags == D4:
            X = (in_bin_medians[8:len(in_bin_medians)-6])
            y = np.log10(red_rms[8:len(red_rms) - 6])
    elif band == 1:
        plt.ylabel(r'$\sigma_%s$' % 'g')
        plt.xlabel(r'$<%s>$' % 'g')
        if mags == D1:
            X = (in_bin_medians[4:len(in_bin_medians)-2])
            y = np.log10(red_rms[4:len(red_rms) - 2])
        elif mags == D2:
            X = (in_bin_medians[8:len(in_bin_medians)-2])
            y = np.log10(red_rms[8:len(red_rms) - 2])
        elif mags == D3:
            X = (in_bin_medians[9:len(in_bin_medians)-1])
            y = np.log10(red_rms[9:len(red_rms)-1])
        elif mags == D4:
            X = (in_bin_medians[0:len(in_bin_medians)-8])
            y = np.log10(red_rms[0:len(red_rms)-8])
    elif band == 2:
        plt.ylabel(r'$\sigma_%s$' % 'r')
        plt.xlabel(r'$<%s>$' % 'r')
        if mags == D1:
            X = (in_bin_medians[9:len(in_bin_medians)-4])
            y = np.log10(red_rms[9:len(red_rms)-4])
        elif mags == D2:
            X = (in_bin_medians[9:len(in_bin_medians)-2])
            y = np.log10(red_rms[9:len(red_rms)-2])
        elif mags == D3:
            X = (in_bin_medians[12:len(in_bin_medians)-2])
            y = np.log10(red_rms[12:len(red_rms)-2])
        elif mags == D4:
            X = (in_bin_medians[17:len(in_bin_medians)-5])
            y = np.log10(red_rms[17:len(red_rms)-5])
    elif band == 3:
        plt.ylabel(r'$\sigma_%s$' % 'i')
        plt.xlabel(r'$<%s>$' % 'i')
        if mags == D1:
            X = (in_bin_medians[12:len(in_bin_medians)-3])
            y = np.log10(red_rms[12:len(red_rms)-3])
        elif mags == D2:
            X = (in_bin_medians[14:len(in_bin_medians)])
            y = np.log10(red_rms[14:len(red_rms)])
        elif mags == D3:
            X = (in_bin_medians[25:len(in_bin_medians)-1])
            y = np.log10(red_rms[25:len(red_rms)-1])
        elif mags == D4:
            X = (in_bin_medians[16:len(in_bin_medians)-5])
            y = np.log10(red_rms[16:len(red_rms)-5])
    elif band == 4:
        plt.ylabel(r'$\sigma_%s$' % 'z')
        plt.xlabel(r'$<%s>$' % 'z')
        if mags == D1:
            X = (in_bin_medians[10:len(in_bin_medians)-6])
            y = np.log10(red_rms[10:len(red_rms)-6])
        elif mags == D2:
            X = (in_bin_medians[15:len(in_bin_medians)-4])
            y = np.log10(red_rms[15:len(red_rms)-4])
        elif mags == D3:
            X = (in_bin_medians[14:len(in_bin_medians)-4])
            y = np.log10(red_rms[14:len(red_rms)-4])
        elif mags == D4:
            X = (in_bin_medians[14:len(in_bin_medians)-4])
            y = np.log10(red_rms[14:len(red_rms)-4])
    
    X = np.asarray(X).reshape(-1,1)
    y = np.asarray(y).reshape(-1,1)
    degree=4
    polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
    polyreg.fit(X,y)

    plt.plot(X,polyreg.predict(X),color="black")
    plt.plot(X,polyreg.predict(X) + np.log10(2),color="red") 
    plt.plot(X,polyreg.predict(X) + np.log10(3),color="orange") 
    plt.plot(X,polyreg.predict(X) + np.log10(5),color="green") 
    
    
    #plt.plot(X,polyreg.predict(X))
    
    '''
    categorized_objects = []
    for i in range(len(medians)):
        y_val = polyreg.predict(np.asarray(medians[i]).reshape(-1,1))
        if np.log10(all_rms[i]) > (y_val + np.log10(5)):
            categorized_objects.append(3)
        elif np.log10(all_rms[i]) > (y_val + np.log10(3)):
            categorized_objects.append(2)
        elif np.log10(all_rms[i]) > (y_val + np.log10(2)):
            categorized_objects.append(1)
        else:
            categorized_objects.append(0)
    '''
    
    
    plt.title("Polynomial regression with degree "+str(degree))
    #degree=61
    #polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
    #polyreg.fit(X,y)
    #plt.plot(X,polyreg.predict(X),color="black")
    #plt.title("Polynomial regression with degree "+str(degree))
    #print("RMSE: " + str(sklearn.metrics.mean_squared_error(y,polyreg.predict(X))))
    #plt.show()
    #return categorized_objects
    #return medians, all_rms
    


# In[7]:


all_categ_vals = []
for i in range(len(all_mags)):
    all_categ_vals.append([[],[],[]])


# In[34]:





#This chunk of code retrieves the rms values by each field and band
categ_1_u = get_categ_vals(0,D1,D1_sigma)
categ_1_g = get_categ_vals(1,D1,D1_sigma)
categ_1_r = get_categ_vals(2,D1,D1_sigma)

categ_2_u = get_categ_vals(0,D2,D2_sigma)
categ_2_g = get_categ_vals(1,D2,D2_sigma)
categ_2_r = get_categ_vals(2,D2,D2_sigma)


categ_3_u = get_categ_vals(0,D3,D3_sigma)
categ_3_g = get_categ_vals(1,D3,D3_sigma)
categ_3_r = get_categ_vals(2,D3,D3_sigma)


categ_4_u = get_categ_vals(0,D4,D4_sigma)
categ_4_g = get_categ_vals(1,D4,D4_sigma)
categ_4_r = get_categ_vals(2,D4,D4_sigma)


# In[36]:


counter_field_1 = -1
counter_field_2 = -1
counter_field_3 = -1
counter_field_4 = -1

for i in range(len(all_mags)):
    if i in fields_1:
        counter_field_1 += 1
        all_categ_vals[i][0] = categ_1_u[counter_field_1]
        all_categ_vals[i][1] = categ_1_g[counter_field_1]
        all_categ_vals[i][2] = categ_1_r[counter_field_1]
        
    elif i in fields_2:
        counter_field_2 += 1
        all_categ_vals[i][0] = categ_2_u[counter_field_2]
        all_categ_vals[i][1] = categ_2_g[counter_field_2]
        all_categ_vals[i][2] = categ_2_r[counter_field_2]
        
    elif i in fields_3:
        counter_field_3 += 1
        all_categ_vals[i][0] = categ_3_u[counter_field_3]
        all_categ_vals[i][1] = categ_3_g[counter_field_3]
        all_categ_vals[i][2] = categ_3_r[counter_field_3]
        
    elif i in fields_4:
        counter_field_4 += 1
        all_categ_vals[i][0] = categ_4_u[counter_field_4]
        all_categ_vals[i][1] = categ_4_g[counter_field_4]
        all_categ_vals[i][2] = categ_4_r[counter_field_4]


# In[6]:


folded_vals = df.values.tolist()


# In[5]:


three_obj_indexes = []
for i in range(len(folded_vals)):
    if (np.mean(folded_vals[i]) == 3):
        three_obj_indexes.append(i)


# In[4]:


print(three_obj_indexes)


# In[ ]:


def plot_unfolded_light_curve(objects):
    """
    Function to create unfolded light curve plots for all objects of interest
    
    Parameters
    ---
    objects: array of indicies of objects of interest in file_list
    
    Returns
    ---
    plots of magnitude vs MJD for all objects
    """

    for star_object in objects:
        plt.figure(figsize = (9, 12))
        plt.xlabel('MJD')
        plt.ylabel('Magnitude')
        plt.gca().invert_yaxis()
        plt.scatter(all_mjd[star_object][0], all_mags[star_object][0], s = 5, c = 'blue', label = 'u')
        plt.scatter(all_mjd[star_object][1], all_mags[star_object][1], s = 5, c = 'green', label = 'g')
        plt.scatter(all_mjd[star_object][2], all_mags[star_object][2], s = 5, c = 'purple', label = 'r')
        plt.scatter(all_mjd[star_object][3], all_mags[star_object][3], s = 5, c = 'gold', label = 'i')
        plt.scatter(all_mjd[star_object][4], all_mags[star_object][4], s = 5, c = 'tab:red', label = 'z')
        plt.legend()
        plt.title(fnames[star_object] + " (" + str(star_object) + ")")
        plt.show()
        


# In[72]:


def data_format (object_number):
    """
    Function that formats the data for one object for Lomb Scargle Periodogram.
    Parameters
    ---
    object_number: int for object ID a.k.a. index in all_mags
    Returns
    ---
    filts: array of bands ex. ['u', 'g', 'g', 'r']
    mags: array of apparent magnitudes that is the same length as filts
    mjd: array of mjd
    dy: array of measurement errors
    """
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


# In[1]:


def lomb_scargle_list (object_list):
    max_over_mean = []
    for cur_object in object_list:
        t, mags, dy, filts = data_format(cur_object)
        periods = np.linspace(0.1, 1.0, 100000)
        LS_multi = periodic.LombScargleMultiband(Nterms_base=1, Nterms_band=0)#initiate structure variable
        LS_multi.fit(t, mags, dy, filts)#input our data
        P_multi = LS_multi.periodogram(periods)#function where input is periods
        max_over_mean.append(np.max(P_multi)/np.mean(P_multi))
    return max_over_mean,object_list


# In[3]:


max_over_mean_threes,object_list = lomb_scargle_list(three_obj_indexes)
print(max_over_mean_threes)
print('separate')
print(max_over_mean)

# In[ ]:


narrowed_objects = []
for i in range(len(max_over_mean)):
    if max_over_mean[i] >= 0.20:
        narrowed_objects.append(object_list[i],max_over_mean[i])
        


