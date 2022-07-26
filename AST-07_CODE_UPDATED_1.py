"""°°°
#### OVERVIEW:
The data required for this project can be found at this link: 
https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/cfhtls/dfspt.html
One must download the DeepVarFull.tar.gz file and drag it/download it into the Jupyter notebook or project to work with it. The data should contain approximately 28000 files to be extracted in which each file represents an astronomical object. There is information about the object within each file, like its magnitude measurements and the filters in which the measurements were taken, etc.. Each of the astronomical objects had several brightness measurements taken in each of the six filters ('U','G','R','I1','I2', and 'Z' are our names for them).
In the code below, you will find that we have combined the fourth and fifth filters (i1 and i2) into a single filter, thus having more magnitude measurements than the other filters.
°°°"""
# |%%--%%| <j8m9d8Ftu9|xurYvUhV0K>

#downloading gatsby
import sys
!{sys.executable} -m pip install gatspy

# |%%--%%| <xurYvUhV0K|vICqjPNjok>

#fields[12]

# |%%--%%| <vICqjPNjok|Z4TitWCBC5>

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









%matplotlib notebook

full_file_list = os.listdir('full') # creates a Python list containing the file paths for every object's time series





# |%%--%%| <Z4TitWCBC5|A8jekuRaNe>
"""°°°
The following code is the basis for all the analysis.
°°°"""
# |%%--%%| <A8jekuRaNe|5SthYwNg8K>

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

    with open('full/'+file_name) as data_file:
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

# |%%--%%| <5SthYwNg8K|MR0Y1bQEx8>

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

    with open('full/'+file_name) as data_file:
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

# |%%--%%| <MR0Y1bQEx8|y6FqGHSHRp>

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

# |%%--%%| <y6FqGHSHRp|vzwqDDEWW7>

num_measurements = get_num_measurements(full_file_list)

# |%%--%%| <vzwqDDEWW7|YrfSItQb2F>

file_list = (np.asarray(full_file_list))[num_measurements>15]

# |%%--%%| <YrfSItQb2F|RyABB4xxu4>

num_measurements_by_filter = get_num_measurements(file_list, by_filter=True)

# |%%--%%| <RyABB4xxu4|NYPHYDYMjY>

#there must be at least 15 measurements per band for an object for us to use it
file_list = file_list[np.all(num_measurements_by_filter>=15, axis=1)]

# |%%--%%| <NYPHYDYMjY|FWhOOb4mF7>



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

    

# |%%--%%| <FWhOOb4mF7|FpBMwaCQaA>

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
    
    

# |%%--%%| <FpBMwaCQaA|5890bXDZQs>

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




# |%%--%%| <5890bXDZQs|Y17F1Sfhne>

#deleting '.mjdmag' from file name of each object
for i, fname in enumerate(fnames):
    fnames[i] = fname.replace(".mjdmag", "")

# |%%--%%| <Y17F1Sfhne|vwOLASyBYL>

for i in range(len(fnames)):
    if fnames[i] == "CFHTLS-VAR-J221527.94-180359.8":
        print(i)

# |%%--%%| <vwOLASyBYL|x9xZo8VqIR>

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

# |%%--%%| <x9xZo8VqIR|ODvpoWN0gM>

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

# |%%--%%| <ODvpoWN0gM|U8G8qcGpAI>

#remove original i1 measurements from 3d list
for i in range(len(all_mags)):
    all_mags[i].remove(all_mags[i][3])
    
for i in range(len(all_mjd)):
    all_mjd[i].remove(all_mjd[i][3])
    
for i in range(len(all_magerrs)):
    all_magerrs[i].remove(all_magerrs[i][3])

# |%%--%%| <U8G8qcGpAI|AO05NgpEPI>

#create a copy of all_mags
all_mags_2 = []
for i in range(len(all_mags)):
    all_mags_2.append([[],[],[],[],[]])


for i in range(len(all_mags)):
    for j in range(len(all_mags[i])):
            all_mags_2[i][j] = copy.deepcopy(all_mags[i][j])

# |%%--%%| <AO05NgpEPI|R3lajYflWl>

all_mags_copy = copy.deepcopy(all_mags)

# |%%--%%| <R3lajYflWl|7AqjUgRKVp>

#5 sigma clip all_mags 
for i in range(len(all_mags)):
    all_mags[i][0] = sigma_clip(all_mags[i][0],sigma=5,maxiters=3,masked=False,copy=False)
    all_mags[i][1] = sigma_clip(all_mags[i][1],sigma=5,maxiters=3,masked=False,copy=False)
    all_mags[i][2] = sigma_clip(all_mags[i][2],sigma=5,maxiters=3,masked=False,copy=False)
    all_mags[i][3] = sigma_clip(all_mags[i][3],sigma=5,maxiters=3,masked=False,copy=False)
    all_mags[i][4] = sigma_clip(all_mags[i][4],sigma=5,maxiters=3,masked=False,copy=False)
    



# |%%--%%| <7AqjUgRKVp|eTHDzmJSMD>

#5 sigma clip all_mags_copy and RETURN MASKED
for i in range(len(all_mags_copy)):
    all_mags_copy[i][0] = sigma_clip(all_mags_copy[i][0],sigma=5,maxiters=3,masked=True,copy=False)
    all_mags_copy[i][1] = sigma_clip(all_mags_copy[i][1],sigma=5,maxiters=3,masked=True,copy=False)
    all_mags_copy[i][2] = sigma_clip(all_mags_copy[i][2],sigma=5,maxiters=3,masked=True,copy=False)
    all_mags_copy[i][3] = sigma_clip(all_mags_copy[i][3],sigma=5,maxiters=3,masked=True,copy=False)
    all_mags_copy[i][4] = sigma_clip(all_mags_copy[i][4],sigma=5,maxiters=3,masked=True,copy=False)
    

# |%%--%%| <eTHDzmJSMD|3VLiLIf5mV>

#deleting corresponding sigma-clipped out data in all_mjd and all_magerrs
for i in range(len(all_mags_copy)):
    for j in range(len(all_mags_copy[i])):
        temp_mask = ma.getmaskarray(all_mags_copy[i][j])
        all_mjd[i][j] = np.delete(all_mjd[i][j], temp_mask)
        all_magerrs[i][j] = np.delete(all_magerrs[i][j], temp_mask)

# |%%--%%| <3VLiLIf5mV|SusjWGq2Ye>

#4.5 sigma clip all_mags_2
for i in range(len(all_mags_2)):
    all_mags_2[i][0] = sigma_clip(all_mags_2[i][0],sigma=4.5,maxiters=3,masked=False,copy=False)
    all_mags_2[i][1] = sigma_clip(all_mags_2[i][1],sigma=4.5,maxiters=3,masked=False,copy=False)
    all_mags_2[i][2] = sigma_clip(all_mags_2[i][2],sigma=4.5,maxiters=3,masked=False,copy=False)
    all_mags_2[i][3] = sigma_clip(all_mags_2[i][3],sigma=4.5,maxiters=3,masked=False,copy=False)
    all_mags_2[i][4] = sigma_clip(all_mags_2[i][4],sigma=4.5,maxiters=3,masked=False,copy=False)

# |%%--%%| <SusjWGq2Ye|jhelxdH6xh>

#separating objects by field
#these lists have indicies of objects
fields_1 = [i for i in range(len(fields)) if fields[i] == "D1"]
fields_2 = [i for i in range(len(fields)) if fields[i] == "D2"]
fields_3 = [i for i in range(len(fields)) if fields[i] == "D3"]
fields_4 = [i for i in range(len(fields)) if fields[i] == "D4"]

# |%%--%%| <jhelxdH6xh|7iUGrXIsC7>

#separating measurements of objects by fields
D1 = [all_mags[i] for i in range(len(all_mags)) if i in fields_1]
D2 = [all_mags[i] for i in range(len(all_mags)) if i in fields_2]
D3 = [all_mags[i] for i in range(len(all_mags)) if i in fields_3]
D4 = [all_mags[i] for i in range(len(all_mags)) if i in fields_4]

# |%%--%%| <7iUGrXIsC7|cQbFZ4Z7a1>

#separating 4.5 sigma clipped measurements of objects by fields
D1_sigma = [all_mags_2[i] for i in range(len(all_mags_2)) if i in fields_1]
D2_sigma = [all_mags_2[i] for i in range(len(all_mags_2)) if i in fields_2]
D3_sigma = [all_mags_2[i] for i in range(len(all_mags_2)) if i in fields_3]
D4_sigma = [all_mags_2[i] for i in range(len(all_mags_2)) if i in fields_4]

# |%%--%%| <cQbFZ4Z7a1|beelPpPGP4>

#sorting file names by field
D1_fnames = [fnames[i] for i in range(len(fnames)) if i in fields_1]
D2_fnames = [fnames[i] for i in range(len(fnames)) if i in fields_2]
D3_fnames = [fnames[i] for i in range(len(fnames)) if i in fields_3]
D4_fnames = [fnames[i] for i in range(len(fnames)) if i in fields_4]

# |%%--%%| <beelPpPGP4|s6PpbgAOWs>

def get_intrinsic_rms_squared(band, magss, mags_sigma):
    """
    Function to plot intrinsic rms for a given field and given filter
    
    Parameters:
    band: 0, 1, 2, 3, 4 (for u, g, r, i, z)
    magss: D1, D2, D3, or D4
    mags_sigma: D1_sigma, D2_sigma, etc.
    
    Returns:
    intrinsic_rms_squared
    median_rms_squared
    plot of intrinsic rms squared vs median mag
    """
    intrinsic_rms_squared = []
    median_rms_squared = []
    mags, rms, in_bin_medians, red_rms = return_mags_rms_medians_red(band, magss, mags_sigma)
    plt.figure(figsize = (9, 12))
    plt.ylabel('Intrinsic RMS Squared')
    if band == 0:
        plt.xlabel(r'$<%s>$' % 'u')
        if magss == D1:
            X = (in_bin_medians[6:len(in_bin_medians)-3])
            y = np.log10(red_rms[6:len(red_rms) - 3])
        elif magss == D2:
            X = (in_bin_medians[6:len(in_bin_medians)-5])
            y = np.log10(red_rms[6:len(red_rms) - 5])
        elif magss == D3:
            X = (in_bin_medians[10:len(in_bin_medians)-6])
            y = np.log10(red_rms[10:len(red_rms) - 6])
        elif magss == D4:
            X = (in_bin_medians[8:len(in_bin_medians)-6])
            y = np.log10(red_rms[8:len(red_rms) - 6])
    elif band == 1:
        plt.xlabel(r'$<%s>$' % 'g')
        if magss == D1:
            X = (in_bin_medians[4:len(in_bin_medians)-2])
            y = np.log10(red_rms[4:len(red_rms) - 2])
        elif magss == D2:
            X = (in_bin_medians[8:len(in_bin_medians)-2])
            y = np.log10(red_rms[8:len(red_rms) - 2])
        elif magss == D3:
            X = (in_bin_medians[9:len(in_bin_medians)-1])
            y = np.log10(red_rms[9:len(red_rms)-1])
        elif magss == D4:
            X = (in_bin_medians[0:len(in_bin_medians)-8])
            y = np.log10(red_rms[0:len(red_rms)-8])
    elif band == 2:
        plt.xlabel(r'$<%s>$' % 'r')
        if magss == D1:
            X = (in_bin_medians[9:len(in_bin_medians)-4])
            y = np.log10(red_rms[9:len(red_rms)-4])
        elif magss == D2:
            X = (in_bin_medians[9:len(in_bin_medians)-2])
            y = np.log10(red_rms[9:len(red_rms)-2])
        elif magss == D3:
            X = (in_bin_medians[12:len(in_bin_medians)-2])
            y = np.log10(red_rms[12:len(red_rms)-2])
        elif magss == D4:
            X = (in_bin_medians[17:len(in_bin_medians)-5])
            y = np.log10(red_rms[17:len(red_rms)-5])
    elif band == 3:
        plt.xlabel(r'$<%s>$' % 'i')
        if magss == D1:
            X = (in_bin_medians[12:len(in_bin_medians)-3])
            y = np.log10(red_rms[12:len(red_rms)-3])
        elif magss == D2:
            X = (in_bin_medians[21:len(in_bin_medians)])
            y = np.log10(red_rms[21:len(red_rms)])
        elif magss == D3:
            X = (in_bin_medians[25:len(in_bin_medians)-1])
            y = np.log10(red_rms[25:len(red_rms)-1])
            plt.ylim(-100, 20)
        elif magss == D4:
            X = (in_bin_medians[16:len(in_bin_medians)-5])
            y = np.log10(red_rms[16:len(red_rms)-5])
    elif band == 4:
        plt.xlabel(r'$<%s>$' % 'z')
        if magss == D1:
            X = (in_bin_medians[10:len(in_bin_medians)-6])
            y = np.log10(red_rms[10:len(red_rms)-6])
        elif magss == D2:
            X = (in_bin_medians[15:len(in_bin_medians)-4])
            y = np.log10(red_rms[15:len(red_rms)-4])
        elif magss == D3:
            X = (in_bin_medians[14:len(in_bin_medians)-4])
            y = np.log10(red_rms[14:len(red_rms)-4])
        elif magss == D4:
            X = (in_bin_medians[14:len(in_bin_medians)-4])
            y = np.log10(red_rms[14:len(red_rms)-4])
    degree=4
    X_1 = X.reshape(-1, 1)
    y_1 = y.reshape(-1, 1)
    polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
    polyreg.fit(X_1,y_1)
    #plt.plot(X,polyreg.predict(X),color="black")
    #plt.title("Polynomial regression with degree "+str(degree))
    #print("RMSE: " + str(sklearn.metrics.mean_squared_error(y,polyreg.predict(X))))
    for i in range(len(mags)):
        mags_arr = mags[i]
        mags_arr = mags_arr.reshape(1, -1)
        intrinsic_rms_squared.append(rms[i]*rms[i]- 10 ** polyreg.predict(mags_arr).item(0) * 10 ** polyreg.predict(mags_arr).item(0))
        median_rms_squared.append(10 ** polyreg.predict(mags_arr).item(0) * 10 ** polyreg.predict(mags_arr).item(0))
    #plt.scatter(mags, intrinsic_rms_squared, s = 0.05)
    #plt.show()
    #return intrinsic_rms_squared, median_rms_squared
    return intrinsic_rms_squared

# |%%--%%| <s6PpbgAOWs|Q27LUeO3ev>

def get_all_mags_format():
    """
    Function to return an empty list with the same dimensions as all_mags
    
    Parameters:
    None
    
    Returns:
    return_list: a list of length 26820 where each element has 5 lists
    """
    return_list = []
    for i in range(len(all_mags)):
        return_list.append([[],[],[],[],[]])
    return return_list

# |%%--%%| <Q27LUeO3ev|klNDwEBwXo>

'''
IMPORTANT NOTE: This code is optional to run. Only run it if you need the 
intrinsic rms squared values for field 1 objects.
----------------------------------------------------------------------------
This code creates a list (field_1_irs) with the same dimensions as all_mags.
that has the intrinsic rms squared value for each object's 5 bands in field 1.
'''
irs_1_u = plot_intrinsic_rms(0,D1,D1_sigma)
irs_1_g = plot_intrinsic_rms(1,D1,D1_sigma)
irs_1_r = plot_intrinsic_rms(2,D1,D1_sigma)
irs_1_i = plot_intrinsic_rms(3,D1,D1_sigma)
irs_1_z = plot_intrinsic_rms(4,D1,D1_sigma)

field_1_irs = get_all_mags_format()

for i in range(len(field_1_irs)):
    field_1_irs[i][0] = irs_1_u[i]
    field_1_irs[i][1] = irs_1_g[i]
    field_1_irs[i][2] = irs_1_r[i]
    field_1_irs[i][3] = irs_1_i[i]
    field_1_irs[i][4] = irs_1_z[i]
    

# |%%--%%| <klNDwEBwXo|mraqXxjiYt>



# |%%--%%| <mraqXxjiYt|OCCYXNytyL>



# |%%--%%| <OCCYXNytyL|Y0evGcZoau>



# |%%--%%| <Y0evGcZoau|JgLlv81AMj>



# |%%--%%| <JgLlv81AMj|dUxo24QkmJ>



# |%%--%%| <dUxo24QkmJ|hFjVv9AIjp>



# |%%--%%| <hFjVv9AIjp|NhVhhKWh1V>



# |%%--%%| <NhVhhKWh1V|0w5Zi5Phx4>



# |%%--%%| <0w5Zi5Phx4|wf2xB0mHkm>



# |%%--%%| <wf2xB0mHkm|kG67OvHCnW>



# |%%--%%| <kG67OvHCnW|j8N6TyiBoL>



# |%%--%%| <j8N6TyiBoL|DSQ0E2s2ln>



# |%%--%%| <DSQ0E2s2ln|4sOE2XQKgt>



# |%%--%%| <4sOE2XQKgt|GBAf7Y7Npw>



# |%%--%%| <GBAf7Y7Npw|C8z4wFNFyZ>



# |%%--%%| <C8z4wFNFyZ|BpB0RgJK8m>



# |%%--%%| <BpB0RgJK8m|nhAttGdTGs>



# |%%--%%| <nhAttGdTGs|DlXWfzEgBD>

#Original Function
def get_error(band, mags, mags_for_median):
    '''
    Function that plots median magnitude vs. decimal log rms for a given filter and given field.
    
    Function that returns a list of rms values per object for a given filter and given band
    
    Parameters
    ---
    band: int for filter (0 = u, 1 = g, ..., 4 = z)
    mags: choose from D1, D2, D3, D4
    mags_for_median: choose from 4.5 sigma clipped mags; choose from D1_sigma, D2_sigma, D3_sigma, D4_sigma
    Returns
    ---
    field_errors: a 1-dimension list of rms values per object in a given field (the field relating to the parameter mags)
    '''
    all_true_mags = []
    median_mags = []
    medians = []
    all_rms = []
    all_rms_two = []
    medians_for_pink = []
    #simplified
    for i in range(len(mags)):
        true_mags = [mags[i][band][j] for j in range(len(mags[i][band]))]
        all_true_mags.append(true_mags)
        temp_mags_for_median = [mags_for_median[i][band][j] for j in range(len(mags_for_median[i][band]))]
        median_mags.append(temp_mags_for_median)
    ####
    for i in range(len(all_true_mags)): #going through 5-sigma clipped mags (i is each object)
        squared_differences = []
        median = np.median(all_true_mags[i]) #median for 1 object’s magnitudes in a specific filter (see: band parameter)
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
    plt.scatter(medians,np.log10(all_rms), s = 0.05)
    plt.scatter(in_bin_medians,np.log10(red_rms), s = 13, c = 'crimson' )#c = for_color, cmap = ‘PuRd’
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
    field_errors = []
    for i in range(len(mags)):
        median_val = np.median(mags[i][band])
        error_polyn_val = polyreg.predict(np.asarray(median_val).reshape(-1,1))
        field_errors.append(error_polyn_val)
    return field_errors

# |%%--%%| <DlXWfzEgBD|oN274EU4O5>

#This chunk of code retrieves the rms values by each field and band
error_1_u = get_error(0,D1,D1_sigma)
error_1_g = get_error(1,D1,D1_sigma)
error_1_r = get_error(2,D1,D1_sigma)
error_1_i = get_error(3,D1,D1_sigma)
error_1_z = get_error(4,D1,D1_sigma)

error_2_u = get_error(0,D2,D2_sigma)
error_2_g = get_error(1,D2,D2_sigma)
error_2_r = get_error(2,D2,D2_sigma)
error_2_i = get_error(3,D2,D2_sigma)
error_2_z = get_error(4,D2,D2_sigma)

error_3_u = get_error(0,D3,D3_sigma)
error_3_g = get_error(1,D3,D3_sigma)
error_3_r = get_error(2,D3,D3_sigma)
error_3_i = get_error(3,D3,D3_sigma)
error_3_z = get_error(4,D3,D3_sigma)

error_4_u = get_error(0,D4,D4_sigma)
error_4_g = get_error(1,D4,D4_sigma)
error_4_r = get_error(2,D4,D4_sigma)
error_4_i = get_error(3,D4,D4_sigma)
error_4_z = get_error(4,D4,D4_sigma)

# |%%--%%| <oN274EU4O5|m2F8AiWEgK>

'''
This code creates a list to hold error values for each object and an object's 5 individual bands.
This list is in the same dimensions as all_mags
'''
all_error_vals = get_all_mags_format()

# |%%--%%| <m2F8AiWEgK|5LsNTbbOQK>

'''
This chunk of code copies the rms values from the several lists (eg. error_1_u) into the single list all_error_vals
by object and band
'''
for i in range(len(all_mags)):
    counter_1 = -1
    counter_2 = -1
    counter_3 = -1
    counter_4 = -1
    if i in fields_1:
        counter_1 += 1
        all_error_vals[i][0] = error_1_u[counter_1]
        all_error_vals[i][1] = error_1_g[counter_1]
        all_error_vals[i][2] = error_1_r[counter_1]
        all_error_vals[i][3] = error_1_i[counter_1]
        all_error_vals[i][4] = error_1_z[counter_1]
    
    elif i in fields_2:
        counter_2 += 1
        all_error_vals[i][0] = error_2_u[counter_2]
        all_error_vals[i][1] = error_2_g[counter_2]
        all_error_vals[i][2] = error_2_r[counter_2]
        all_error_vals[i][3] = error_2_i[counter_2]
        all_error_vals[i][4] = error_2_z[counter_2]
        
    elif i in fields_3:
        counter_3 += 1
        all_error_vals[i][0] = error_3_u[counter_3]
        all_error_vals[i][1] = error_3_g[counter_3]
        all_error_vals[i][2] = error_3_r[counter_3]
        all_error_vals[i][3] = error_3_i[counter_3]
        all_error_vals[i][4] = error_3_z[counter_3]
        
    
    elif i in fields_4:
        counter_4 += 1
        all_error_vals[i][0] = error_4_u[counter_4]
        all_error_vals[i][1] = error_4_g[counter_4]
        all_error_vals[i][2] = error_4_r[counter_4]
        all_error_vals[i][3] = error_4_i[counter_4]
        all_error_vals[i][4] = error_4_z[counter_2]
        
    else:
        print("Not working!")

for i in range(len(all_error_vals)):
    for j in range(len(all_error_vals[i])):
        all_error_vals[i][j] = all_error_vals[i][j][0][0]

# |%%--%%| <5LsNTbbOQK|U9OKvPnS1g>

#This cell of code creates a new list (all_percentages) to hold the rms error percentages per object and band
all_percentages = get_all_mags_format()
for i in range(len(all_error_vals)):
    for j in range(len(all_error_vals[i])):
        neg_val = -1 * all_error_vals[i][j]
        error_val = pow(10,neg_val)
        range_val = abs(max(all_mags[i][j]) - min(all_mags[i][j]))
        percentage = error_val / range_val
        all_percentages[i][j] = percentage

# |%%--%%| <U9OKvPnS1g|HlGKFs9And>
"""°°°
all_percentages holds the rms error percentage for each object's 5 bands
(eg. all_percentages[0][0] is the rms error percentage for the zeroth object in 
the u band).
°°°"""
# |%%--%%| <HlGKFs9And|68IkWwSfKB>
"""°°°
The following code runs the lomb-scargle periodogram and creates folded light curves.
°°°"""
# |%%--%%| <68IkWwSfKB|TQipABwTaB>

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

# |%%--%%| <TQipABwTaB|nQfzTomlwr>

t, mags, dy, filts = data_format(15615)

# |%%--%%| <nQfzTomlwr|5Htq76235G>

#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib as mpl
#mpl.style.use('ggplot')
#mpl.rc('axes', color_cycle=["#4C72B0", "#55A868", "#C44E52","#8172B2", "#CCB974"]) colors for each band
#from gatspy.periodic import LombScargleMultiband
#from gatspy import datasets, periodic

# Choose a Sesar 2010 object to base our fits on
#lcid = 1019544
#rrlyrae = datasets.RRLyraeGenerated(lcid, random_state=0)

# Generate data in a 6-month observing season
#Nobs = 60
#rng = np.random.RandomState(0)

#nights = np.arange(180)
#rng.shuffle(nights)
#nights = nights[:Nobs]

#t = 57000 + nights + 0.05 * rng.randn(Nobs)
#dy = 0.06 + 0.01 * rng.randn(Nobs)
#mags = np.array([rrlyrae.generated(band, t, err=dy, corrected=False)for band in 'ugriz'])
"""
All the above are not necessary for the NGVS_legacy work, becuase you do have realistic data.
I generate mock datesets just to test to code and show an example result.
"""

"""
Here start the fitting process. 
t, mags, dy, filts are all N-d arraies of your observation data. They mean observation time in MJD, apperent magnitudes, errors, and filter list, respectively
If you do not have error list, set dy to be all ones.
filts could be anything in the format of np.array(['u','u','g', ..., 'z']), as long as its elements are filter names and its length equal to the mags array
"""
#example of appropriate format
#filts = np.take(list('ugriz'), np.arange(Nobs), mode='wrap') 
# 
#mags = mags[np.arange(Nobs) % 5, np.arange(Nobs)]
#masks = [(filts == band) for band in 'ugriz']#separate ugriz to 5 sublists

#below is the necessary code
periods = np.linspace(0.1, 1.0, 100000) # This defines the search range of your period, you can specify it at your will. These are in days.


#2 different ways of fitting
##model = periodic.NaiveMultiband(BaseModel=periodic.LombScargleFast) 
# specify the method to be the naive multiband LS, which means you fit data in each band separately, and get a score list for each band.
# serves as a good first try on your NGVS_legacy data
#above is the fastest way. x axis is period, y axis is score. 5d array output. real variable should have same peak for all bands
##model.fit(t, mags, dy, filts) 
##P = model.scores(periods) 
# This is the fitting score list you want. 
# It is a 5xN array, where N is number of periods tried, P[0] is the fit socres of periods with your u band data. And so on for P[1] for g, P[2] for i, ...
# Each element in P[i] correspond to a period in the array periods you input. The closer to 1, the better.


#all bands done together
LS_multi = periodic.LombScargleMultiband(Nterms_base=1, Nterms_band=0)#initiate structure variable
LS_multi.fit(t, mags, dy, filts)#input our data
P_multi = LS_multi.periodogram(periods)#function where input is periods
# A non-naive way of multiband fitting. This time all data from all bands are fitted simultaneously, means you do not get scores for each band separately.
# P_multi will be a 1-d array, has equal length to your input periods. The maximum value in P_multi will be the best fit, and its corresponding period will be the best period.

#do both!!

"""
From here are visualization of the results.

fig = plt.figure(figsize=(10, 4))
gs = plt.GridSpec(5, 2, left=0.07, right=0.95, bottom=0.15,
                  wspace=0.1, hspace=0.6)
ax = [fig.add_subplot(gs[:, 0]),
      fig.add_subplot(gs[:-2, 1]),
      fig.add_subplot(gs[-2:, 1])]

for band, mask in zip('ugriz', masks):
    ax[0].errorbar((t[mask] / rrlyrae.period) % 1, mags[mask], dy[mask],
                   fmt='.', label=band)
ax[0].set_ylim(18, 14.5)
ax[0].legend(loc='upper left', fontsize=12, ncol=3)
ax[0].set_title('Folded Data, 1 band per night (P={0:.3f} days)'
                ''.format(rrlyrae.period), fontsize=12)
ax[0].set_xlabel('phase')
ax[0].set_ylabel('magnitude')

for i, band in enumerate('ugriz'):
    offset = 4 - i
    ax[1].plot(periods, P[band] + offset, lw=1)
    ax[1].text(0.89, 1 + offset, band, fontsize=10, ha='right', va='top')
ax[1].set_title('Standard Periodogram in Each Band', fontsize=12)
ax[1].yaxis.set_major_formatter(plt.NullFormatter())
ax[1].xaxis.set_major_formatter(plt.NullFormatter())
ax[1].set_ylabel('power + offset')


ax[2].plot(periods, P_multi, lw=1, color='gray')

ax[2].set_title('Multiband Periodogram', fontsize=12)
ax[2].set_yticks([0, 0.5, 1.0])
ax[2].set_ylim(0, 1.0)
ax[2].yaxis.set_major_formatter(plt.NullFormatter())
ax[2].set_xlabel('Period (days)')
ax[2].set_ylabel('power')
"""

# |%%--%%| <5Htq76235G|58W2zRIyQA>

plt.ion() #turns figure display on
plt.figure()
plt.scatter(periods, P_multi, s = 0.05)

# |%%--%%| <58W2zRIyQA|olKkVqDrNi>

best_period = max(P_multi)
for i in range(len(P_multi)):
    if P_multi[i] == best_period:
        index = i
print(periods[index])

# |%%--%%| <olKkVqDrNi|DsKMjAkPu2>

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
        #for star_object in objects:
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
        plt.title(fnames[obj] + " (" + str(obj) + ")")
        pdf.savefig()
        plt.close()
    #0.6 seconds ideal step

# |%%--%%| <DsKMjAkPu2|0n0nDXtek0>

folded_light_curve(503, periods[index], fnames[503] + "million.pdf")

# |%%--%%| <0n0nDXtek0|Hzy5TjCSCA>

def not_pdf_folded_light_curve(obj, period):
    #rr lyrae has variation in g band of 0.6 or 0.5 in unfolded light curve
    plt.ion()
    #with PdfPages(pdf_name) as pdf:
        #for star_object in objects:
    plt.figure(figsize = (9, 12))
    plt.xlabel('Phase')
    plt.ylabel('Magnitude')
    plt.gca().invert_yaxis()
    peak_mag = min(all_mags[obj][0])
    u_phase = [(i%period)/period for i in all_mjd[obj][0]]
    g_phase = [(i%period)/period for i in all_mjd[obj][1]]
    r_phase = [(i%period)/period for i in all_mjd[obj][2]]
    i_phase = [(i%period)/period for i in all_mjd[obj][3]]
    z_phase = [(i%period)/period for i in all_mjd[obj][4]]
    for i in range(len(all_mags[obj][0])):
        if all_mags[obj][0][i] == peak_mag:
            peak_index = i
    u_phase_change = [u_phase[i]- u_phase[peak_index] if u_phase[i]- u_phase[peak_index] >= 0 else u_phase[i]- u_phase[peak_index] + 1 for i in range(len(u_phase))]
    g_phase_change = [g_phase[i]- u_phase[peak_index] if g_phase[i]- u_phase[peak_index] >= 0 else g_phase[i]- u_phase[peak_index] + 1 for i in range(len(g_phase))]
    r_phase_change = [r_phase[i]- u_phase[peak_index] if r_phase[i]- u_phase[peak_index] >= 0 else r_phase[i]- u_phase[peak_index] + 1 for i in range(len(r_phase))]
    i_phase_change = [i_phase[i]- u_phase[peak_index] if i_phase[i]- u_phase[peak_index] >= 0 else i_phase[i]- u_phase[peak_index] + 1 for i in range(len(i_phase))]
    z_phase_change = [z_phase[i]- u_phase[peak_index] if z_phase[i]- u_phase[peak_index] >= 0 else z_phase[i]- u_phase[peak_index] + 1 for i in range(len(z_phase))]
        #f*p^2/length of observations, f is about 0.1, f = phase error for individual cycle
    plt.scatter(u_phase_change, all_mags[obj][0], s = 5, c = 'blue', label = 'u')
    plt.scatter(g_phase_change, all_mags[obj][1], s = 5, c = 'green', label = 'g')
    plt.scatter(r_phase_change, all_mags[obj][2], s = 5, c = 'purple', label = 'r')
    plt.scatter(i_phase_change, all_mags[obj][3], s = 5, c = 'gold', label = 'i')
    plt.scatter(z_phase_change, all_mags[obj][4], s = 5, c = 'tab:red', label = 'z')
    plt.legend()
    plt.title(fnames[obj] + " (" + str(obj) + ")")
        #pdf.savefig()
        #plt.close()
    #0.6 seconds ideal step
    # i/p_true + n where n is an integer all reciprocated is the beat frequency

# |%%--%%| <Hzy5TjCSCA|eczhRN91i9>

def plot_rms_mags (band, mags, mags_for_median):
    “”"
    Function that plots median magnitude vs. decimal log rms for a given filter and given field.
    Parameters
    ---
    band: int for filter (0 = u, 1 = g, ..., 4 = z)
    mags: choose from D1, D2, D3, D4
    mags_for_median: choose from 4.5 sigma clipped mags; choose from D1_sigma, D2_sigma, D3_sigma, D4_sigma
    Returns
    ---
    Plot of median magnitudes vs. log rms for a specific band and specific field
    “”"
    all_true_mags = []
    median_mags = []
    medians = []
    all_rms = []
    all_rms_two = []
    medians_for_pink = []
    “”"
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
    “”"
    #simplified
    for i in range(len(mags)):
        true_mags = [mags[i][band][j] for j in range(len(mags[i][band]))]
        all_true_mags.append(true_mags)
        temp_mags_for_median = [mags_for_median[i][band][j] for j in range(len(mags_for_median[i][band]))]
        median_mags.append(temp_mags_for_median)
    ####
    for i in range(len(all_true_mags)): #going through 5-sigma clipped mags (i is each object)
        squared_differences = []
        median = np.median(all_true_mags[i]) #median for 1 object’s magnitudes in a specific filter (see: band parameter)
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
    plt.scatter(medians,np.log10(all_rms), s = 0.05)
    plt.scatter(in_bin_medians,np.log10(red_rms), s = 13, c = “crimson”)#c = for_color, cmap = ‘PuRd’
    if band == 0:
        plt.ylabel(r’$\sigma_%s$' % ‘u’)
        plt.xlabel(r’$<%s>$' % ‘u’)
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
        plt.ylabel(r’$\sigma_%s$' % ‘g’)
        plt.xlabel(r’$<%s>$' % ‘g’)
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
        plt.ylabel(r’$\sigma_%s$' % ‘r’)
        plt.xlabel(r’$<%s>$' % ‘r’)
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
        plt.ylabel(r’$\sigma_%s$' % ‘i’)
        plt.xlabel(r’$<%s>$' % ‘i’)
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
        plt.ylabel(r’$\sigma_%s$' % ‘z’)
        plt.xlabel(r’$<%s>$' % ‘z’)
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
    plt.plot(X,polyreg.predict(X),color=“black”)
    plt.title(“Polynomial regression with degree “+str(degree))
    #degree=61
    #polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
    #polyreg.fit(X,y)
    #plt.plot(X,polyreg.predict(X),color=“black”)
    #plt.title(“Polynomial regression with degree “+str(degree))
    #print(“RMSE: ” + str(sklearn.metrics.mean_squared_error(y,polyreg.predict(X))))
    plt.show()
    #return medians, all_rms
                   
                   
                   
    #Optional: get medians
    medians_obj_band = all_mags.copy()
    for i in range(len(all_mags)):
        for j in range(len(all_mags[i])):
            medians_obj_band[i][j] = np.median(all_mags[i][j])
            
    Y_X = None
    error_measurements = medians_obj_band.copy()
    for i in range(len(medians_obj_band)):
        for j in range(len(medians_obj_band[i])):
            Y_X = polyreg.predict(np.asarray(medians_obj_band[i][j]).reshape(-1,1))
            error_measurements[i][j] = Y_X
            
            
    error_measurements_list = all_mags.copy()
    for i in range(len(error_measurements_list)):
        for j in range(5):
            error_measurements_list[i][j] = error_measurements[i][j][0][0]

        
    return error_measurements_list
                

# |%%--%%| <eczhRN91i9|q2RFrYtegl>

not_pdf_folded_light_curve(15615, periods[index])

# |%%--%%| <q2RFrYtegl|ah3ROPPora>

plt.figure()
plt.hist(all_magerrs[503][1], bins = 50)
plt.show()

# |%%--%%| <ah3ROPPora|6ZBCPXM89L>
"""°°°
The following code is used to categorize objects into marginal, intermediate, extreme, or non variables.
°°°"""
# |%%--%%| <6ZBCPXM89L|l09VJOzCPB>

def plot_rms_mags (band, mags, mags_for_median):
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
    plt.scatter(medians,np.log10(all_rms), s = 0.05)
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
    plt.title("Polynomial regression with degree "+str(degree))
    #degree=61
    #polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
    #polyreg.fit(X,y)
    #plt.plot(X,polyreg.predict(X),color="black")
    #plt.title("Polynomial regression with degree "+str(degree))
    #print("RMSE: " + str(sklearn.metrics.mean_squared_error(y,polyreg.predict(X))))
    plt.show()
    #return medians, all_rms
    
    #Optional: get medians
    medians_obj_band = all_mags.copy()
    for i in range(len(all_mags)):
        for j in range(len(all_mags[i])):
            medians_obj_band[i][j] = np.median(all_mags[i][j])
            
    Y_X = None
    error_measurements = medians_obj_band.copy()
    for i in range(len(medians_obj_band)):
        for j in range(len(medians_obj_band[i])):
            Y_X = polyreg.predict(np.asarray(medians_obj_band[i][j]).reshape(-1,1))
            error_measurements[i][j] = Y_X
            
            
    error_measurements_list = all_mags.copy()
    for i in range(len(error_measurements_list)):
        for j in range(5):
            error_measurements_list[i][j] = error_measurements[i][j][0][0]

        
    return error_measurements_list

# |%%--%%| <l09VJOzCPB|vpm7NB6T17>

D1_u_medians, D1_u_rms = plot_rms_mags(0, D1, D1_sigma)

# |%%--%%| <vpm7NB6T17|RrJ91UI5VX>

D1_g_medians, D1_g_rms = plot_rms_mags(1, D1, D1_sigma)

# |%%--%%| <RrJ91UI5VX|EIj9EysVwz>

D1_r_medians, D1_r_rms = plot_rms_mags(2, D1, D1_sigma)

# |%%--%%| <EIj9EysVwz|zC53kN5qfH>

D1_i_medians, D1_i_rms = plot_rms_mags(3, D1, D1_sigma)

# |%%--%%| <zC53kN5qfH|NZP7Yc2BGD>

D1_z_medians, D1_z_rms = plot_rms_mags(4, D1, D1_sigma)

# |%%--%%| <NZP7Yc2BGD|EywvuBfrLR>

D2_u_medians, D2_u_rms = plot_rms_mags(0, D2, D2_sigma)

# |%%--%%| <EywvuBfrLR|lwUbZRlSgS>

D2_g_medians, D2_g_rms = plot_rms_mags(1, D2, D2_sigma)

# |%%--%%| <lwUbZRlSgS|OQC1UORmUP>

D2_r_medians, D2_r_rms = plot_rms_mags(2, D2, D2_sigma)

# |%%--%%| <OQC1UORmUP|onj4LcUrS5>

D2_i_medians, D2_i_rms = plot_rms_mags(3, D2, D2_sigma)

# |%%--%%| <onj4LcUrS5|xgZbaYAa7n>

D2_z_medians, D2_z_rms = plot_rms_mags(4, D2, D2_sigma)

# |%%--%%| <xgZbaYAa7n|W2fbCB44KE>

D3_u_medians, D3_u_rms = plot_rms_mags(0, D3, D3_sigma)

# |%%--%%| <W2fbCB44KE|rEzBhDs3he>

D3_g_medians, D3_g_rms = plot_rms_mags(1, D3, D3_sigma)

# |%%--%%| <rEzBhDs3he|LGditWQuZ8>

D3_r_medians, D3_r_rms = plot_rms_mags(2, D3, D3_sigma)

# |%%--%%| <LGditWQuZ8|GhN1Vq4IP5>

D3_i_medians, D3_i_rms = plot_rms_mags(3, D3, D3_sigma)

# |%%--%%| <GhN1Vq4IP5|n8fW3I7B0A>

D3_z_medians, D3_z_rms = plot_rms_mags(4, D3, D3_sigma)

# |%%--%%| <n8fW3I7B0A|qdmdVckoQG>

D4_u_medians, D4_u_rms = plot_rms_mags(0, D4, D4_sigma)

# |%%--%%| <qdmdVckoQG|qQ7qAwtriy>

D4_g_medians, D4_g_rms = plot_rms_mags(1, D4, D4_sigma)

# |%%--%%| <qQ7qAwtriy|q01r0UaX5w>

D4_r_medians, D4_r_rms = plot_rms_mags(2, D4, D4_sigma)

# |%%--%%| <q01r0UaX5w|mflPR7r3Rq>

D4_i_medians, D4_i_rms = plot_rms_mags(3, D4, D4_sigma)

# |%%--%%| <mflPR7r3Rq|T3QvUp32jC>

D4_z_medians, D4_z_rms = plot_rms_mags(4, D4, D4_sigma)

# |%%--%%| <T3QvUp32jC|s3SG3b4o3f>

D1_u_rms_squared = np.square(D1_u_rms)
D1_g_rms_squared = np.square(D1_g_rms)
D1_r_rms_squared = np.square(D1_r_rms)
D1_i_rms_squared = np.square(D1_i_rms)
D1_z_rms_squared = np.square(D1_z_rms)

D2_u_rms_squared = np.square(D2_u_rms)
D2_g_rms_squared = np.square(D2_g_rms)
D2_r_rms_squared = np.square(D2_r_rms)
D2_i_rms_squared = np.square(D2_i_rms)
D2_z_rms_squared = np.square(D2_z_rms)

D3_u_rms_squared = np.square(D3_u_rms)
D3_g_rms_squared = np.square(D3_g_rms)
D3_r_rms_squared = np.square(D3_r_rms)
D3_i_rms_squared = np.square(D3_i_rms)
D3_z_rms_squared = np.square(D3_z_rms)

D4_u_rms_squared = np.square(D4_u_rms)
D4_g_rms_squared = np.square(D4_g_rms)
D4_r_rms_squared = np.square(D4_r_rms)
D4_i_rms_squared = np.square(D4_i_rms)
D4_z_rms_squared = np.square(D4_z_rms)

# |%%--%%| <s3SG3b4o3f|JrG81oYP40>

def not_plot_rms_mags (band, mags, mags_for_median, err):
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
    #plt.figure(figsize = (9, 12))
    #plt.scatter(medians,np.log10(all_rms), s = 0.05)
    #plt.scatter(in_bin_medians,np.log10(red_rms), s = 13, c = "crimson")#c = for_color, cmap = 'PuRd'
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
    #plt.plot(X,polyreg.predict(X),color="black")
    #plt.title("Polynomial regression with degree "+str(degree))
    #degree=61
    #polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
    #polyreg.fit(X,y)
    #plt.plot(X,polyreg.predict(X),color="black")
    #plt.title("Polynomial regression with degree "+str(degree))
    #print("RMSE: " + str(sklearn.metrics.mean_squared_error(y,polyreg.predict(X))))
    val = np.array(err)
    val_input = val.reshape(1, -1)
    val_output = polyreg.predict(val_input)
    #plt.show()
    return 10**val_output

# |%%--%%| <JrG81oYP40|zLcUEAVzm1>

another_stuff = not_plot_rms_mags(1, D4, D4_sigma, 24)

# |%%--%%| <zLcUEAVzm1|s9cIRRsDAf>

another_stuff

# |%%--%%| <s9cIRRsDAf|jHVXjQ4ofI>

10**another_stuff

# |%%--%%| <jHVXjQ4ofI|WobtlFpVWX>

#Does what plot_rms_mags does EXCEPT the x axis is decimal log rms and the y axis is median magnitude
#Used for finding the median magnitude given a measurement error
def plot_rms_mags_inverted (band, mags, mags_for_median, err):
    """
    Function that plots median magnitude vs. decimal log rms for a given filter and given field.
    Parameters
    ---
    band: int for filter (0 = u, 1 = g, ..., 4 = z)
    mags: choose from D1, D2, D3, D4
    mags_for_median: choose from D1_sigma, D2_sigma, D3_sigma, D4_sigma
    err: desired measurement error (NOT IN DECIMAL LOG)
    Returns
    ---
    var_output: median magnitude for specified err
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
    #plt.figure(figsize = (9, 12))
    #plt.scatter(np.log10(all_rms),medians, s = 0.05)
    #plt.scatter(np.log10(red_rms),in_bin_medians, s = 13, c = "crimson")#c = for_color, cmap = 'PuRd'
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
    polyreg.fit(y,X)
    #FORNOWplt.plot(y,polyreg.predict(y),color="black")
    plt.title("Polynomial regression with degree "+str(degree))
    parameter = np.log10(err)
    var_input = parameter.reshape(1, -1)
    var_output = polyreg.predict(var_input)
    #degree=61
    #polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
    #polyreg.fit(X,y)
    #plt.plot(X,polyreg.predict(X),color="black")
    #plt.title("Polynomial regression with degree "+str(degree))
    #print("RMSE: " + str(sklearn.metrics.mean_squared_error(y,polyreg.predict(X))))
    #FORNOWplt.show()
    return var_output

# |%%--%%| <WobtlFpVWX|LXcfxzsTOz>

stuff = plot_rms_mags_inverted(4, D4, D4_sigma, 0.2)

# |%%--%%| <LXcfxzsTOz|zle6NsYKAC>

print(stuff)

# |%%--%%| <zle6NsYKAC|SbqOcEe3EJ>

def plot_rms_mags_together(band, mags, mags_for_median, field_name, color_val = "black"):
    """
    Function that plots median magnitude vs. rms for a given filter.
    Parameters
    ---
    band: int for filter (0 = u, 1 = g, ..., 4 = z)
    mags: all_mags
    Returns
    ---
    Plot of median magnitudes vs. rms
    """
    
    all_true_mags = []
    median_mags = []
    medians = []
    all_rms = []
    all_rms_two = []
    medians_for_pink = []
    
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
    #plt.figure(figsize = (9, 12))
    #plt.scatter(medians,all_rms, s = 0.05, c = "steelblue")
    plt.scatter(in_bin_medians,red_rms, s = 13, c = "steelblue")#c = for_color, cmap = 'PuRd'
    if band == 0:
        plt.ylabel(r'$\sigma_%s$' % 'u')
        plt.xlabel(r'$<%s>$' % 'u')
        if mags == D1:
            X = (in_bin_medians[6:len(in_bin_medians)-3])
            y = (red_rms[6:len(red_rms) - 3])
        elif mags == D2:
            X = (in_bin_medians[6:len(in_bin_medians)-5])
            y = (red_rms[6:len(red_rms) - 5])
        elif mags == D3:
            X = (in_bin_medians[10:len(in_bin_medians)-6])
            y = (red_rms[10:len(red_rms) - 6])
        elif mags == D4:
            X = (in_bin_medians[8:len(in_bin_medians)-6])
            y = (red_rms[8:len(red_rms) - 6])
    elif band == 1:
        plt.ylabel(r'$\sigma_%s$' % 'g')
        plt.xlabel(r'$<%s>$' % 'g')
        if mags == D1:
            X = (in_bin_medians[4:len(in_bin_medians)-2])
            y = (red_rms[4:len(red_rms) - 2])
        elif mags == D2:
            X = (in_bin_medians[8:len(in_bin_medians)-2])
            y = (red_rms[8:len(red_rms) - 2])
        elif mags == D3:
            X = (in_bin_medians[9:len(in_bin_medians)-1])
            y = (red_rms[9:len(red_rms)-1])
        elif mags == D4:
            X = (in_bin_medians[0:len(in_bin_medians)-8])
            y = (red_rms[0:len(red_rms)-8])
    elif band == 2:
        plt.ylabel(r'$\sigma_%s$' % 'r')
        plt.xlabel(r'$<%s>$' % 'r')
        if mags == D1:
            X = (in_bin_medians[9:len(in_bin_medians)-4])
            y = (red_rms[9:len(red_rms)-4])
        elif mags == D2:
            X = (in_bin_medians[9:len(in_bin_medians)-2])
            y = (red_rms[9:len(red_rms)-2])
        elif mags == D3:
            X = (in_bin_medians[12:len(in_bin_medians)-2])
            y = (red_rms[12:len(red_rms)-2])
        elif mags == D4:
            X = (in_bin_medians[17:len(in_bin_medians)-5])
            y = (red_rms[17:len(red_rms)-5])
    elif band == 3:
        plt.ylabel(r'$\sigma_%s$' % 'i')
        plt.xlabel(r'$<%s>$' % 'i')
        if mags == D1:
            X = (in_bin_medians[12:len(in_bin_medians)-3])
            y = (red_rms[12:len(red_rms)-3])
        elif mags == D2:
            X = (in_bin_medians[14:len(in_bin_medians)])
            y = (red_rms[14:len(red_rms)])
        elif mags == D3:
            X = (in_bin_medians[25:len(in_bin_medians)-1])
            y = (red_rms[25:len(red_rms)-1])
        elif mags == D4:
            X = (in_bin_medians[16:len(in_bin_medians)-5])
            y = (red_rms[16:len(red_rms)-5])
    elif band == 4:
        plt.ylabel(r'$\sigma_%s$' % 'z')
        plt.xlabel(r'$<%s>$' % 'z')
        if mags == D1:
            X = (in_bin_medians[10:len(in_bin_medians)-6])
            y = (red_rms[10:len(red_rms)-6])
        elif mags == D2:
            X = (in_bin_medians[15:len(in_bin_medians)-4])
            y = (red_rms[15:len(red_rms)-4])
        elif mags == D3:
            X = (in_bin_medians[14:len(in_bin_medians)-4])
            y = (red_rms[14:len(red_rms)-4])
        elif mags == D4:
            X = (in_bin_medians[14:len(in_bin_medians)-4])
            y = (red_rms[14:len(red_rms)-4])
    
    X = np.asarray(X).reshape(-1,1)
    y = np.asarray(y).reshape(-1,1)
    degree=4
    polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
    polyreg.fit(X,y)
    plt.plot(X,polyreg.predict(X),color=color_val, label = field_name)
    plt.title("Polynomial regression with degree "+str(degree))
    #degree=61
    #polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
    #polyreg.fit(X,y)
    #plt.plot(X,polyreg.predict(X),color="black")
    #plt.title("Polynomial regression with degree "+str(degree))
    #print("RMSE: " + str(sklearn.metrics.mean_squared_error(y,polyreg.predict(X))))
    plt.show()

# |%%--%%| <SbqOcEe3EJ|XvwK9eN4lV>

#Plots median mag vs rms with ALL 4 fields for a specified band
all_D = [D1, D2, D3, D4]
all_D_sigma = [D1_sigma, D2_sigma, D3_sigma, D4_sigma]
color_vals = ["red", "gold", "purple", "green"]
field_names = ["D1", "D2", "D3", "D4"]
plt.figure(figsize = (9, 12))
for i in range(4):
    plot_rms_mags_together(0, all_D[i], all_D_sigma[i], field_names[i], color_vals[i]) #change the first parameter to change the band
plt.legend(loc = "lower right")

# |%%--%%| <XvwK9eN4lV|lTTRwI4JIY>

#NEEDED TO RUN plot_intrinsic_rms
def return_mags_rms_medians_red (band, mags, mags_for_median):
    """
    Function that returns information for plot_rms_mags.
    Parameters
    ---
    band: int for filter (0 = u, 1 = g, ..., 4 = z)
    mags: D1, D2, D3, or D4
    mags_for_median: D1_sigma, D2_sigma, etc.
    
    Returns
    ---
    median mag, rms of mags, medians sorted by 0.2 mag bins, and rms sorted by 0.2 mag bins
    """
    all_true_mags = []
    median_mags = []
    medians = []
    all_rms = []
    all_rms_two = []
    medians_for_pink = []
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
    #print(all_rms_two)
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
    return medians, all_rms, in_bin_medians, red_rms







# |%%--%%| <lTTRwI4JIY|cewyqCiGPm>

def plot_intrinsic_rms(band, magss, mags_sigma):
    """
    Function to plot intrinsic rms for a given field and given filter
    
    Parameters:
    band: 0, 1, 2, 3, 4 (for u, g, r, i, z)
    magss: D1, D2, D3, or D4
    mags_sigma: D1_sigma, D2_sigma, etc.
    
    Returns:
    intrinsic_rms_squared
    median_rms_squared
    plot of intrinsic rms squared vs median mag
    """
    intrinsic_rms_squared = []
    median_rms_squared = []
    mags, rms, in_bin_medians, red_rms = return_mags_rms_medians_red(band, magss, mags_sigma)
    plt.figure(figsize = (9, 12))
    plt.ylabel('Intrinsic RMS Squared')
    if band == 0:
        plt.xlabel(r'$<%s>$' % 'u')
        if magss == D1:
            X = (in_bin_medians[6:len(in_bin_medians)-3])
            y = np.log10(red_rms[6:len(red_rms) - 3])
        elif magss == D2:
            X = (in_bin_medians[6:len(in_bin_medians)-5])
            y = np.log10(red_rms[6:len(red_rms) - 5])
        elif magss == D3:
            X = (in_bin_medians[10:len(in_bin_medians)-6])
            y = np.log10(red_rms[10:len(red_rms) - 6])
        elif magss == D4:
            X = (in_bin_medians[8:len(in_bin_medians)-6])
            y = np.log10(red_rms[8:len(red_rms) - 6])
    elif band == 1:
        plt.xlabel(r'$<%s>$' % 'g')
        if magss == D1:
            X = (in_bin_medians[4:len(in_bin_medians)-2])
            y = np.log10(red_rms[4:len(red_rms) - 2])
        elif magss == D2:
            X = (in_bin_medians[8:len(in_bin_medians)-2])
            y = np.log10(red_rms[8:len(red_rms) - 2])
        elif magss == D3:
            X = (in_bin_medians[9:len(in_bin_medians)-1])
            y = np.log10(red_rms[9:len(red_rms)-1])
        elif magss == D4:
            X = (in_bin_medians[0:len(in_bin_medians)-8])
            y = np.log10(red_rms[0:len(red_rms)-8])
    elif band == 2:
        plt.xlabel(r'$<%s>$' % 'r')
        if magss == D1:
            X = (in_bin_medians[9:len(in_bin_medians)-4])
            y = np.log10(red_rms[9:len(red_rms)-4])
        elif magss == D2:
            X = (in_bin_medians[9:len(in_bin_medians)-2])
            y = np.log10(red_rms[9:len(red_rms)-2])
        elif magss == D3:
            X = (in_bin_medians[12:len(in_bin_medians)-2])
            y = np.log10(red_rms[12:len(red_rms)-2])
        elif magss == D4:
            X = (in_bin_medians[17:len(in_bin_medians)-5])
            y = np.log10(red_rms[17:len(red_rms)-5])
    elif band == 3:
        plt.xlabel(r'$<%s>$' % 'i')
        if magss == D1:
            X = (in_bin_medians[12:len(in_bin_medians)-3])
            y = np.log10(red_rms[12:len(red_rms)-3])
        elif magss == D2:
            X = (in_bin_medians[21:len(in_bin_medians)])
            y = np.log10(red_rms[21:len(red_rms)])
        elif magss == D3:
            X = (in_bin_medians[25:len(in_bin_medians)-1])
            y = np.log10(red_rms[25:len(red_rms)-1])
            plt.ylim(-100, 20)
        elif magss == D4:
            X = (in_bin_medians[16:len(in_bin_medians)-5])
            y = np.log10(red_rms[16:len(red_rms)-5])
    elif band == 4:
        plt.xlabel(r'$<%s>$' % 'z')
        if magss == D1:
            X = (in_bin_medians[10:len(in_bin_medians)-6])
            y = np.log10(red_rms[10:len(red_rms)-6])
        elif magss == D2:
            X = (in_bin_medians[15:len(in_bin_medians)-4])
            y = np.log10(red_rms[15:len(red_rms)-4])
        elif magss == D3:
            X = (in_bin_medians[14:len(in_bin_medians)-4])
            y = np.log10(red_rms[14:len(red_rms)-4])
        elif magss == D4:
            X = (in_bin_medians[14:len(in_bin_medians)-4])
            y = np.log10(red_rms[14:len(red_rms)-4])
    degree=4
    X_1 = X.reshape(-1, 1)
    y_1 = y.reshape(-1, 1)
    polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
    polyreg.fit(X_1,y_1)
    #plt.plot(X,polyreg.predict(X),color="black")
    #plt.title("Polynomial regression with degree "+str(degree))
    #print("RMSE: " + str(sklearn.metrics.mean_squared_error(y,polyreg.predict(X))))
    for i in range(len(mags)):
        mags_arr = mags[i]
        mags_arr = mags_arr.reshape(1, -1)
        intrinsic_rms_squared.append(rms[i]*rms[i]- 10 ** polyreg.predict(mags_arr).item(0) * 10 ** polyreg.predict(mags_arr).item(0))
        median_rms_squared.append(10 ** polyreg.predict(mags_arr).item(0) * 10 ** polyreg.predict(mags_arr).item(0))
    plt.scatter(mags, intrinsic_rms_squared, s = 0.05)
    plt.show()
    return intrinsic_rms_squared, median_rms_squared

# |%%--%%| <cewyqCiGPm|9jIdu4WJSk>

D1_u_intrinsic_rms_squared, D1_u_median_rms_squared = plot_intrinsic_rms(0, D1, D1_sigma)

# |%%--%%| <9jIdu4WJSk|nHBkrLdcG2>

D1_g_intrinsic_rms_squared, D1_g_median_rms_squared = plot_intrinsic_rms(1, D1, D1_sigma)

# |%%--%%| <nHBkrLdcG2|EiCwnNlmD0>

D1_r_intrinsic_rms_squared, D1_r_median_rms_squared = plot_intrinsic_rms(2, D1, D1_sigma)

# |%%--%%| <EiCwnNlmD0|QMdwj7topN>

D1_i_intrinsic_rms_squared, D1_i_median_rms_squared = plot_intrinsic_rms(3, D1, D1_sigma)

# |%%--%%| <QMdwj7topN|2BgLwBY7nx>

D1_z_intrinsic_rms_squared, D1_z_median_rms_squared = plot_intrinsic_rms(4, D1, D1_sigma)

# |%%--%%| <2BgLwBY7nx|zd5STihZ46>

D2_u_intrinsic_rms_squared, D2_u_median_rms_squared = plot_intrinsic_rms(0, D2, D2_sigma)

# |%%--%%| <zd5STihZ46|9GmaUsXOmj>

D2_g_intrinsic_rms_squared, D2_g_median_rms_squared = plot_intrinsic_rms(1, D2, D2_sigma)

# |%%--%%| <9GmaUsXOmj|hyPd5P0gur>

D2_r_intrinsic_rms_squared, D2_r_median_rms_squared = plot_intrinsic_rms(2, D2, D2_sigma)

# |%%--%%| <hyPd5P0gur|tEbp0iHPOm>

D2_i_intrinsic_rms_squared, D2_i_median_rms_squared = plot_intrinsic_rms(3, D2, D2_sigma)

# |%%--%%| <tEbp0iHPOm|zSPZrYEY8O>

D2_z_intrinsic_rms_squared, D2_z_median_rms_squared = plot_intrinsic_rms(4, D2, D2_sigma)

# |%%--%%| <zSPZrYEY8O|lSHUXPvvN7>

D3_u_intrinsic_rms_squared, D3_u_median_rms_squared = plot_intrinsic_rms(0, D3, D3_sigma)

# |%%--%%| <lSHUXPvvN7|MkQtMzeMRW>

D3_g_intrinsic_rms_squared, D3_g_median_rms_squared = plot_intrinsic_rms(1, D3, D3_sigma)

# |%%--%%| <MkQtMzeMRW|y0ifVLH8nM>

D3_r_intrinsic_rms_squared, D3_r_median_rms_squared = plot_intrinsic_rms(2, D3, D3_sigma)

# |%%--%%| <y0ifVLH8nM|4phufg8lLN>

D3_i_intrinsic_rms_squared, D3_i_median_rms_squared = plot_intrinsic_rms(3, D3, D3_sigma)

# |%%--%%| <4phufg8lLN|pcn1yNfiKw>

D3_z_intrinsic_rms_squared, D3_z_median_rms_squared = plot_intrinsic_rms(4, D3, D3_sigma)

# |%%--%%| <pcn1yNfiKw|NlZWlnBSUU>

D4_u_intrinsic_rms_squared, D4_u_median_rms_squared = plot_intrinsic_rms(0, D4, D4_sigma)

# |%%--%%| <NlZWlnBSUU|v9gXvci1QL>

D4_g_intrinsic_rms_squared, D4_g_median_rms_squared = plot_intrinsic_rms(1, D4, D4_sigma)

# |%%--%%| <v9gXvci1QL|BYgTqSsHZz>

D4_r_intrinsic_rms_squared, D4_r_median_rms_squared = plot_intrinsic_rms(2, D4, D4_sigma)

# |%%--%%| <BYgTqSsHZz|XAT3ij0hG7>

D4_i_intrinsic_rms_squared, D4_i_median_rms_squared = plot_intrinsic_rms(3, D4, D4_sigma)

# |%%--%%| <XAT3ij0hG7|WCnmfukmVP>

D4_z_intrinsic_rms_squared, D4_z_median_rms_squared = plot_intrinsic_rms(4, D4, D4_sigma)

# |%%--%%| <WCnmfukmVP|xwjINqa72S>

#create boolean to make sure median is above the cutoff due to the tail on the left hand side of rms vs median mag
D1_z_medians_boolean = [np.median(D1[i][4]) > 16.5 for i in range(len(D1))]
D2_z_medians_boolean = [np.median(D2[i][4]) > 16.5 for i in range(len(D2))]
D3_z_medians_boolean = [np.median(D3[i][4]) > 16.5 for i in range(len(D3))]
D4_z_medians_boolean = [np.median(D4[i][4]) > 16.5 for i in range(len(D4))]

D1_i_medians_boolean = [np.median(D1[i][3]) > 18 for i in range(len(D1))]
D2_i_medians_boolean = [np.median(D2[i][3]) > 18 for i in range(len(D2))]
D3_i_medians_boolean = [np.median(D3[i][3]) > 18 for i in range(len(D3))]
D4_i_medians_boolean = [np.median(D4[i][3]) > 18 for i in range(len(D4))]

D1_r_medians_boolean = [np.median(D1[i][2]) > 17.5 for i in range(len(D1))]
D2_r_medians_boolean = [np.median(D2[i][2]) > 17.5 for i in range(len(D2))]
D3_r_medians_boolean = [np.median(D3[i][2]) > 17.5 for i in range(len(D3))]
D4_r_medians_boolean = [np.median(D4[i][2]) > 17.5 for i in range(len(D4))]

# |%%--%%| <xwjINqa72S|xSIqKAdkAT>

#reciprocals for weights for weighted summed intrinsic rms squared 
D1_u_reciprocal = np.reciprocal(D1_u_median_rms_squared)
D1_g_reciprocal = np.reciprocal(D1_g_median_rms_squared)
D1_r_reciprocal = np.reciprocal(D1_r_median_rms_squared)
D1_i_reciprocal = np.reciprocal(D1_i_median_rms_squared)
D1_z_reciprocal = np.reciprocal(D1_z_median_rms_squared)
D2_u_reciprocal = np.reciprocal(D2_u_median_rms_squared)
D2_g_reciprocal = np.reciprocal(D2_g_median_rms_squared)
D2_r_reciprocal = np.reciprocal(D2_r_median_rms_squared)
D2_i_reciprocal = np.reciprocal(D2_i_median_rms_squared)
D2_z_reciprocal = np.reciprocal(D2_z_median_rms_squared)
D3_u_reciprocal = np.reciprocal(D3_u_median_rms_squared)
D3_g_reciprocal = np.reciprocal(D3_g_median_rms_squared)
D3_r_reciprocal = np.reciprocal(D3_r_median_rms_squared)
D3_i_reciprocal = np.reciprocal(D3_i_median_rms_squared)
D3_z_reciprocal = np.reciprocal(D3_z_median_rms_squared)
D4_u_reciprocal = np.reciprocal(D4_u_median_rms_squared)
D4_g_reciprocal = np.reciprocal(D4_g_median_rms_squared)
D4_r_reciprocal = np.reciprocal(D4_r_median_rms_squared)
D4_i_reciprocal = np.reciprocal(D4_i_median_rms_squared)
D4_z_reciprocal = np.reciprocal(D4_z_median_rms_squared)

D1_sum_intrinsic_rms = []
D2_sum_intrinsic_rms = []
D3_sum_intrinsic_rms = []
D4_sum_intrinsic_rms = []

for i in range(len(D1_u_intrinsic_rms_squared)):
    numerator = D1_u_reciprocal[i] * D1_u_intrinsic_rms_squared[i] + D1_g_reciprocal[i] * D1_g_intrinsic_rms_squared[i] + D1_r_reciprocal[i] * D1_r_intrinsic_rms_squared[i] + D1_i_reciprocal[i] * D1_i_intrinsic_rms_squared[i] + D1_z_reciprocal[i] * D1_z_intrinsic_rms_squared[i]
    denominator = D1_u_reciprocal[i] + D1_g_reciprocal[i] + D1_r_reciprocal[i] + D1_i_reciprocal[i] + D1_z_reciprocal[i]
    D1_sum_intrinsic_rms.append(numerator/denominator)

for i in range(len(D2_u_intrinsic_rms_squared)):
    numerator = D2_u_reciprocal[i] * D2_u_intrinsic_rms_squared[i] + D2_g_reciprocal[i] * D2_g_intrinsic_rms_squared[i] + D2_r_reciprocal[i] * D2_r_intrinsic_rms_squared[i] + D2_i_reciprocal[i] * D2_i_intrinsic_rms_squared[i] + D2_z_reciprocal[i] * D2_z_intrinsic_rms_squared[i]
    denominator = D2_u_reciprocal[i] + D2_g_reciprocal[i] + D2_r_reciprocal[i] + D2_i_reciprocal[i] + D2_z_reciprocal[i]
    D2_sum_intrinsic_rms.append(numerator/denominator)

for i in range(len(D3_u_intrinsic_rms_squared)):
    numerator = D3_u_reciprocal[i] * D3_u_intrinsic_rms_squared[i] + D3_g_reciprocal[i] * D3_g_intrinsic_rms_squared[i] + D3_r_reciprocal[i] * D3_r_intrinsic_rms_squared[i] + D3_i_reciprocal[i] * D3_i_intrinsic_rms_squared[i] + D3_z_reciprocal[i] * D3_z_intrinsic_rms_squared[i]
    denominator = D3_u_reciprocal[i] + D3_g_reciprocal[i] + D3_r_reciprocal[i] + D3_i_reciprocal[i] + D3_z_reciprocal[i]
    D3_sum_intrinsic_rms.append(numerator/denominator)

for i in range(len(D4_u_intrinsic_rms_squared)):
    numerator = D4_u_reciprocal[i] * D4_u_intrinsic_rms_squared[i] + D4_g_reciprocal[i] * D4_g_intrinsic_rms_squared[i] + D4_r_reciprocal[i] * D4_r_intrinsic_rms_squared[i] + D4_i_reciprocal[i] * D4_i_intrinsic_rms_squared[i] + D4_z_reciprocal[i] * D4_z_intrinsic_rms_squared[i]
    denominator = D4_u_reciprocal[i] + D4_g_reciprocal[i] + D4_r_reciprocal[i] + D4_i_reciprocal[i] + D4_z_reciprocal[i]
    D4_sum_intrinsic_rms.append(numerator/denominator)

# |%%--%%| <xSIqKAdkAT|Z6jJUqQBgl>

#to find index of object in file_list ERROR
"""
D1_indicies = []
D2_indicies = []
D3_indicies = []
D4_indicies = []
for i in range(len(D1)):
    D1_indicies.append(i)
for i in range(len(D2)):
    D2_indicies.append(i + len(D1))
for i in range(len(D3)):
    D3_indicies.append(i + len(D1) + len(D2))
for i in range(len(D4)):
    D4_indicies.append(i + len(D1) + len(D2) + len(D3))
"""

# |%%--%%| <Z6jJUqQBgl|0zMPGS3LNN>

#testing to see if the objects satisfy the median mag requirements/cutoff
D1_indicies_boolean = [D1_z_medians_boolean[i] == True and D1_r_medians_boolean[i] == True and D1_i_medians_boolean[i] == True for i in range(len(D1_z_medians_boolean))]
D2_indicies_boolean = [D2_z_medians_boolean[i] == True and D2_r_medians_boolean[i] == True and D2_i_medians_boolean[i] == True for i in range(len(D2_z_medians_boolean))]
D3_indicies_boolean = [D3_z_medians_boolean[i] == True and D3_r_medians_boolean[i] == True and D3_i_medians_boolean[i] == True for i in range(len(D3_z_medians_boolean))]
D4_indicies_boolean = [D4_z_medians_boolean[i] == True and D4_r_medians_boolean[i] == True and D4_i_medians_boolean[i] == True for i in range(len(D4_z_medians_boolean))]

# |%%--%%| <0zMPGS3LNN|kASnT5R3w3>

#changes true to false and vice versa
D1_indicies_boolean = [not elem for elem in D1_indicies_boolean]
D2_indicies_boolean = [not elem for elem in D2_indicies_boolean]
D3_indicies_boolean = [not elem for elem in D3_indicies_boolean]
D4_indicies_boolean = [not elem for elem in D4_indicies_boolean]

# |%%--%%| <kASnT5R3w3|izOJy2dfNr>

#do not use objects that don't satisfy median mag cutoff
D1 = np.delete(fields_1, D1_indicies_boolean)
D2 = np.delete(fields_2, D2_indicies_boolean)
D3 = np.delete(fields_3, D3_indicies_boolean)
D4 = np.delete(fields_4, D4_indicies_boolean)

# |%%--%%| <izOJy2dfNr|PKHg4YiWyU>

#delete file names of those objects that are not used
D1_fnames = np.delete(D1_fnames, D1_indicies_boolean)
D2_fnames = np.delete(D2_fnames, D2_indicies_boolean)
D3_fnames = np.delete(D3_fnames, D3_indicies_boolean)
D4_fnames = np.delete(D4_fnames, D4_indicies_boolean)

# |%%--%%| <PKHg4YiWyU|9R3Ct6hcbF>

#delete RMS of objects that are not used
D1_u_rms = np.delete(D1_u_rms, D1_indicies_boolean)
D1_g_rms = np.delete(D1_g_rms, D1_indicies_boolean)
D1_r_rms = np.delete(D1_r_rms, D1_indicies_boolean)
D1_i_rms = np.delete(D1_i_rms, D1_indicies_boolean)
D1_z_rms = np.delete(D1_z_rms, D1_indicies_boolean)

D2_u_rms = np.delete(D2_u_rms, D2_indicies_boolean)
D2_g_rms = np.delete(D2_g_rms, D2_indicies_boolean)
D2_r_rms = np.delete(D2_r_rms, D2_indicies_boolean)
D2_i_rms = np.delete(D2_i_rms, D2_indicies_boolean)
D2_z_rms = np.delete(D2_z_rms, D2_indicies_boolean)

D3_u_rms = np.delete(D3_u_rms, D3_indicies_boolean)
D3_g_rms = np.delete(D3_g_rms, D3_indicies_boolean)
D3_r_rms = np.delete(D3_r_rms, D3_indicies_boolean)
D3_i_rms = np.delete(D3_i_rms, D3_indicies_boolean)
D3_z_rms = np.delete(D3_z_rms, D3_indicies_boolean)

D4_u_rms = np.delete(D4_u_rms, D4_indicies_boolean)
D4_g_rms = np.delete(D4_g_rms, D4_indicies_boolean)
D4_r_rms = np.delete(D4_r_rms, D4_indicies_boolean)
D4_i_rms = np.delete(D4_i_rms, D4_indicies_boolean)
D4_z_rms = np.delete(D4_z_rms, D4_indicies_boolean)

# |%%--%%| <9R3Ct6hcbF|7bdyUCIm1f>

#delete sum intrinsic rms of those objects that are not used
D1_sum_intrinsic_rms = [D1_sum_intrinsic_rms[i] for i in range(len(D1_sum_intrinsic_rms)) if D1_z_medians_boolean[i] == True and D1_i_medians_boolean[i] == True and D1_r_medians_boolean[i] == True]
D2_sum_intrinsic_rms = [D2_sum_intrinsic_rms[i] for i in range(len(D2_sum_intrinsic_rms)) if D2_z_medians_boolean[i] == True and D2_i_medians_boolean[i] == True and D2_r_medians_boolean[i] == True]
D3_sum_intrinsic_rms = [D3_sum_intrinsic_rms[i] for i in range(len(D3_sum_intrinsic_rms)) if D3_z_medians_boolean[i] == True and D3_i_medians_boolean[i] == True and D3_r_medians_boolean[i] == True]
D4_sum_intrinsic_rms = [D4_sum_intrinsic_rms[i] for i in range(len(D4_sum_intrinsic_rms)) if D4_z_medians_boolean[i] == True and D4_i_medians_boolean[i] == True and D4_r_medians_boolean[i] == True]

# |%%--%%| <7bdyUCIm1f|GYt3Xo7UMf>

#finding median g mags for each field
D1_g_median_mags = [np.median(D1_sigma[i][1]) for i in range(len(D1_sigma)) if D1_z_medians_boolean[i] == True and D1_i_medians_boolean[i] == True and D1_r_medians_boolean[i] == True]
D2_g_median_mags = [np.median(D2_sigma[i][1]) for i in range(len(D2_sigma)) if D2_z_medians_boolean[i] == True and D2_i_medians_boolean[i] == True and D2_r_medians_boolean[i] == True]
D3_g_median_mags = [np.median(D3_sigma[i][1]) for i in range(len(D3_sigma)) if D3_z_medians_boolean[i] == True and D3_i_medians_boolean[i] == True and D3_r_medians_boolean[i] == True]
D4_g_median_mags = [np.median(D4_sigma[i][1]) for i in range(len(D4_sigma)) if D4_z_medians_boolean[i] == True and D4_i_medians_boolean[i] == True and D4_r_medians_boolean[i] == True]

# |%%--%%| <GYt3Xo7UMf|bfiYTg2ena>

D1_g_bin_medians = []
D1_g_bin_rms = []
D2_g_bin_medians = []
D2_g_bin_rms = []
D3_g_bin_medians = []
D3_g_bin_rms = []
D4_g_bin_medians = []
D4_g_bin_rms = []

# |%%--%%| <bfiYTg2ena|3RJoRJpBfH>

#creating mag bins
D1_temp_bin_1 = [D1_g_median_mags[i] for i in range(len(D1_g_median_mags)) if 17 <= D1_g_median_mags[i] < 18 ]
D1_temp_bin_1a = [D1_sum_intrinsic_rms[i] for i in range(len(D1_g_median_mags)) if 17 <= D1_g_median_mags[i] < 18]
D1_temp_bin_2 = [D1_g_median_mags[i] for i in range(len(D1_g_median_mags)) if 18 <= D1_g_median_mags[i] < 19]
D1_temp_bin_2a = [D1_sum_intrinsic_rms[i] for i in range(len(D1_g_median_mags)) if 18 <= D1_g_median_mags[i] < 19]
D1_temp_bin_3 = [D1_g_median_mags[i] for i in range(len(D1_g_median_mags)) if 19 <= D1_g_median_mags[i] < 19.5]
D1_temp_bin_3a = [D1_sum_intrinsic_rms[i] for i in range(len(D1_g_median_mags)) if 19 <= D1_g_median_mags[i] < 19.5]
D1_temp_bin_4 = [D1_g_median_mags[i] for i in range(len(D1_g_median_mags)) if 19.5 <= D1_g_median_mags[i] < 20]
D1_temp_bin_4a = [D1_sum_intrinsic_rms[i] for i in range(len(D1_g_median_mags)) if 19.5 <= D1_g_median_mags[i] < 20]

D2_temp_bin_1 = [D2_g_median_mags[i] for i in range(len(D2_g_median_mags)) if 17 <= D2_g_median_mags[i] < 18 ]
D2_temp_bin_1a = [D2_sum_intrinsic_rms[i] for i in range(len(D2_g_median_mags)) if 17 <= D2_g_median_mags[i] < 18]
D2_temp_bin_2 = [D2_g_median_mags[i] for i in range(len(D2_g_median_mags)) if 18 <= D2_g_median_mags[i] < 19]
D2_temp_bin_2a = [D2_sum_intrinsic_rms[i] for i in range(len(D2_g_median_mags)) if 18 <= D2_g_median_mags[i] < 19]
D2_temp_bin_3 = [D2_g_median_mags[i] for i in range(len(D2_g_median_mags)) if 19 <= D2_g_median_mags[i] < 19.5]
D2_temp_bin_3a = [D2_sum_intrinsic_rms[i] for i in range(len(D2_g_median_mags)) if 19 <= D2_g_median_mags[i] < 19.5]
D2_temp_bin_4 = [D2_g_median_mags[i] for i in range(len(D2_g_median_mags)) if 19.5 <= D2_g_median_mags[i] < 20]
D2_temp_bin_4a = [D2_sum_intrinsic_rms[i] for i in range(len(D2_g_median_mags)) if 19.5 <= D2_g_median_mags[i] < 20]

D3_temp_bin_1 = [D3_g_median_mags[i] for i in range(len(D3_g_median_mags)) if 17 <= D3_g_median_mags[i] < 18 ]
D3_temp_bin_1a = [D3_sum_intrinsic_rms[i] for i in range(len(D3_g_median_mags)) if 17 <= D3_g_median_mags[i] < 18]
D3_temp_bin_2 = [D3_g_median_mags[i] for i in range(len(D3_g_median_mags)) if 18 <= D3_g_median_mags[i] < 19]
D3_temp_bin_2a = [D3_sum_intrinsic_rms[i] for i in range(len(D3_g_median_mags)) if 18 <= D3_g_median_mags[i] < 19]
D3_temp_bin_3 = [D3_g_median_mags[i] for i in range(len(D3_g_median_mags)) if 19 <= D3_g_median_mags[i] < 19.5]
D3_temp_bin_3a = [D3_sum_intrinsic_rms[i] for i in range(len(D3_g_median_mags)) if 19 <= D3_g_median_mags[i] < 19.5]
D3_temp_bin_4 = [D3_g_median_mags[i] for i in range(len(D3_g_median_mags)) if 19.5 <= D3_g_median_mags[i] < 20]
D3_temp_bin_4a = [D3_sum_intrinsic_rms[i] for i in range(len(D3_g_median_mags)) if 19.5 <= D3_g_median_mags[i] < 20]

D4_temp_bin_1 = [D4_g_median_mags[i] for i in range(len(D4_g_median_mags)) if 17 <= D4_g_median_mags[i] < 18 ]
D4_temp_bin_1a = [D4_sum_intrinsic_rms[i] for i in range(len(D4_g_median_mags)) if 17 <= D4_g_median_mags[i] < 18]
D4_temp_bin_2 = [D4_g_median_mags[i] for i in range(len(D4_g_median_mags)) if 18 <= D4_g_median_mags[i] < 19]
D4_temp_bin_2a = [D4_sum_intrinsic_rms[i] for i in range(len(D4_g_median_mags)) if 18 <= D4_g_median_mags[i] < 19]
D4_temp_bin_3 = [D4_g_median_mags[i] for i in range(len(D4_g_median_mags)) if 19 <= D4_g_median_mags[i] < 19.5]
D4_temp_bin_3a = [D4_sum_intrinsic_rms[i] for i in range(len(D4_g_median_mags)) if 19 <= D4_g_median_mags[i] < 19.5]
D4_temp_bin_4 = [D4_g_median_mags[i] for i in range(len(D4_g_median_mags)) if 19.5 <= D4_g_median_mags[i] < 20]
D4_temp_bin_4a = [D4_sum_intrinsic_rms[i] for i in range(len(D4_g_median_mags)) if 19.5 <= D4_g_median_mags[i] < 20]

# |%%--%%| <3RJoRJpBfH|A5AxgOQRxy>

#appending bins to a list
D1_g_bin_medians.append(D1_temp_bin_1)
D1_g_bin_medians.append(D1_temp_bin_2)
D1_g_bin_medians.append(D1_temp_bin_3)
D1_g_bin_medians.append(D1_temp_bin_4)
D1_g_bin_rms.append(D1_temp_bin_1a)
D1_g_bin_rms.append(D1_temp_bin_2a)
D1_g_bin_rms.append(D1_temp_bin_3a)
D1_g_bin_rms.append(D1_temp_bin_4a)

D2_g_bin_medians.append(D2_temp_bin_1)
D2_g_bin_medians.append(D2_temp_bin_2)
D2_g_bin_medians.append(D2_temp_bin_3)
D2_g_bin_medians.append(D2_temp_bin_4)
D2_g_bin_rms.append(D2_temp_bin_1a)
D2_g_bin_rms.append(D2_temp_bin_2a)
D2_g_bin_rms.append(D2_temp_bin_3a)
D2_g_bin_rms.append(D2_temp_bin_4a)

D3_g_bin_medians.append(D3_temp_bin_1)
D3_g_bin_medians.append(D3_temp_bin_2)
D3_g_bin_medians.append(D3_temp_bin_3)
D3_g_bin_medians.append(D3_temp_bin_4)
D3_g_bin_rms.append(D3_temp_bin_1a)
D3_g_bin_rms.append(D3_temp_bin_2a)
D3_g_bin_rms.append(D3_temp_bin_3a)
D3_g_bin_rms.append(D3_temp_bin_4a)

D4_g_bin_medians.append(D4_temp_bin_1)
D4_g_bin_medians.append(D4_temp_bin_2)
D4_g_bin_medians.append(D4_temp_bin_3)
D4_g_bin_medians.append(D4_temp_bin_4)
D4_g_bin_rms.append(D4_temp_bin_1a)
D4_g_bin_rms.append(D4_temp_bin_2a)
D4_g_bin_rms.append(D4_temp_bin_3a)
D4_g_bin_rms.append(D4_temp_bin_4a)

# |%%--%%| <A5AxgOQRxy|SyrsGh9Wev>

#appending 0.2 mag bins to list
bounds = np.arange(20, 24.6, 0.2)
for i in range(len(bounds) -1):
    temp_bin_medians = []
    temp_bin_rms = []
    for j in range(len(D1_g_median_mags)):
        if bounds[i] <= D1_g_median_mags[j] < bounds[i+1]:
            temp_bin_medians.append(D1_g_median_mags[j])
            temp_bin_rms.append(D1_sum_intrinsic_rms[j])
    D1_g_bin_medians.append(temp_bin_medians)
    D1_g_bin_rms.append(temp_bin_rms)
    
for i in range(len(bounds) -1):
    temp_bin_medians = []
    temp_bin_rms = []
    for j in range(len(D2_g_median_mags)):
        if bounds[i] <= D2_g_median_mags[j] < bounds[i+1]:
            temp_bin_medians.append(D2_g_median_mags[j])
            temp_bin_rms.append(D2_sum_intrinsic_rms[j])
    D2_g_bin_medians.append(temp_bin_medians)
    D2_g_bin_rms.append(temp_bin_rms)
    
for i in range(len(bounds) -1):
    temp_bin_medians = []
    temp_bin_rms = []
    for j in range(len(D3_g_median_mags)):
        if bounds[i] <= D3_g_median_mags[j] < bounds[i+1]:
            temp_bin_medians.append(D3_g_median_mags[j])
            temp_bin_rms.append(D3_sum_intrinsic_rms[j])
    D3_g_bin_medians.append(temp_bin_medians)
    D3_g_bin_rms.append(temp_bin_rms)
    
for i in range(len(bounds) -1):
    temp_bin_medians = []
    temp_bin_rms = []
    for j in range(len(D4_g_median_mags)):
        if bounds[i] <= D4_g_median_mags[j] < bounds[i+1]:
            temp_bin_medians.append(D4_g_median_mags[j])
            temp_bin_rms.append(D4_sum_intrinsic_rms[j])
    D4_g_bin_medians.append(temp_bin_medians)
    D4_g_bin_rms.append(temp_bin_rms)

# |%%--%%| <SyrsGh9Wev|GtTMe0yHzr>

#copying medians just in case deleting goes wrong
D1_g_bin_medians_copy = copy.deepcopy(D1_g_bin_medians)
D2_g_bin_medians_copy = copy.deepcopy(D2_g_bin_medians)
D3_g_bin_medians_copy = copy.deepcopy(D3_g_bin_medians)
D4_g_bin_medians_copy = copy.deepcopy(D4_g_bin_medians)

# |%%--%%| <GtTMe0yHzr|pu6JQYUlVu>

#some rms are 0 which causes issues, so deleting those
for bin_index in range(len(D1_g_bin_rms)):
    D1_g_bin_rms[bin_index] = [i for i in D1_g_bin_rms[bin_index] if i < 0]

for bin_index in range(len(D2_g_bin_rms)):
    D2_g_bin_rms[bin_index] = [i for i in D2_g_bin_rms[bin_index] if i < 0]
    
for bin_index in range(len(D3_g_bin_rms)):
    D3_g_bin_rms[bin_index] = [i for i in D3_g_bin_rms[bin_index] if i < 0]
    
for bin_index in range(len(D4_g_bin_rms)):
    D4_g_bin_rms[bin_index] = [i for i in D4_g_bin_rms[bin_index] if i < 0]

# |%%--%%| <pu6JQYUlVu|nu2lrMGCY1>

#creating boolean to test for empty arrays
D1_boolean = [len(D1_g_bin_rms[i]) == 0 for i in range(len(D1_g_bin_rms))]
D2_boolean = [len(D2_g_bin_rms[i]) == 0 for i in range(len(D2_g_bin_rms))]
D3_boolean = [len(D3_g_bin_rms[i]) == 0 for i in range(len(D3_g_bin_rms))]
D4_boolean = [len(D4_g_bin_rms[i]) == 0 for i in range(len(D4_g_bin_rms))]

# |%%--%%| <nu2lrMGCY1|aKcRHFjPlF>

#deleting empty arrays
D1_g_bin_rms = np.delete(D1_g_bin_rms, D1_boolean)
D1_g_bin_medians = np.delete(D1_g_bin_medians, D1_boolean)

D2_g_bin_rms = np.delete(D2_g_bin_rms, D2_boolean)
D2_g_bin_medians = np.delete(D2_g_bin_medians, D2_boolean)

D3_g_bin_rms = np.delete(D3_g_bin_rms, D3_boolean)
D3_g_bin_medians = np.delete(D3_g_bin_medians, D3_boolean)

D4_g_bin_rms = np.delete(D4_g_bin_rms, D4_boolean)
D4_g_bin_medians = np.delete(D4_g_bin_medians, D4_boolean)

# |%%--%%| <aKcRHFjPlF|1fXCztwkis>

#finding 99.98th percentile rms for each field
D1_x_percentile_vals = []
D1_y_percentile_vals = []
D2_x_percentile_vals = []
D2_y_percentile_vals = []
D3_x_percentile_vals = []
D3_y_percentile_vals = []
D4_x_percentile_vals = []
D4_y_percentile_vals = []
for i in range(len(D1_g_bin_rms)):
    D1_y_percentile_vals.append(np.percentile(np.absolute(D1_g_bin_rms[i]), 99.98))
    D1_x_percentile_vals.append(np.median(D1_g_bin_medians[i]))

for i in range(len(D2_g_bin_rms)):
    D2_y_percentile_vals.append(np.percentile(np.absolute(D2_g_bin_rms[i]), 99.98))
    D2_x_percentile_vals.append(np.median(D2_g_bin_medians[i]))
    
for i in range(len(D3_g_bin_rms)):
    D3_y_percentile_vals.append(np.percentile(np.absolute(D3_g_bin_rms[i]), 99.98))
    D3_x_percentile_vals.append(np.median(D3_g_bin_medians[i]))
    
for i in range(len(D4_g_bin_rms)):
    D4_y_percentile_vals.append(np.percentile(np.absolute(D4_g_bin_rms[i]), 99.98))
    D4_x_percentile_vals.append(np.median(D4_g_bin_medians[i]))

# |%%--%%| <1fXCztwkis|y8ohIh1kJ5>

#for i in range(len(y_percentile_vals)):
    #y_percentile_vals[i] = y_percentile_vals[i] * -1

# |%%--%%| <y8ohIh1kJ5|Q67HT0IMac>

#FOR D1
X = D1_x_percentile_vals#[:len(x_percentile_vals) -1]
y = D1_y_percentile_vals#[:len(y_percentile_vals)-1]
#print(X)
#print(y)
X = np.asarray(X).reshape(-1, 1)
y = np.asarray(y).reshape(-1, 1)
degree=80
polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
polyreg.fit(X,y)
plt.figure()
plt.scatter(D1_g_median_mags, D1_sum_intrinsic_rms, s = 0.01)
plt.xlabel('Median <g> Magnitudes')
plt.ylabel('Weighted Summed Intrinsic RMS Squared')
plt.xlim(17, 24.5)
plt.ylim(-0.1, 1)
plt.plot(X,polyreg.predict(X), color="black")
plt.axhline(y = 0.0025, color = "crimson")#cutoff line
plt.axhline(y = 0.01, color = "blue")#cutoff line
plt.scatter(D1_x_percentile_vals, D1_y_percentile_vals, s = 5, c = 'red')
plt.title("Polynomial regression with degree "+str(degree))
plt.show()

# |%%--%%| <Q67HT0IMac|MRbKx8zjr6>

#FOR D2
X = D2_x_percentile_vals#[:len(x_percentile_vals) -1]
y = D2_y_percentile_vals#[:len(y_percentile_vals)-1]
#print(X)
#print(y)
X = np.asarray(X).reshape(-1, 1)
y = np.asarray(y).reshape(-1, 1)
degree=71
polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
polyreg.fit(X,y)
plt.figure()
plt.scatter(D2_g_median_mags, D2_sum_intrinsic_rms, s = 0.01)
plt.xlabel('Median <g> Magnitudes')
plt.ylabel('Weighted Summed Intrinsic RMS Squared')
plt.xlim(17, 24.5)
plt.ylim(-0.5, 1)
plt.plot(X,polyreg.predict(X), color="black")
plt.axhline(y = 0.0025, color = "crimson")
plt.axhline(y = 0.01, color = "blue")
plt.scatter(D2_x_percentile_vals, D2_y_percentile_vals, s = 5, c = 'red')
plt.title("Polynomial regression with degree "+str(degree))
plt.show()

# |%%--%%| <MRbKx8zjr6|tmHlxZxdRA>

#FOR D3
X = D3_x_percentile_vals#[:len(x_percentile_vals) -1]
y = D3_y_percentile_vals#[:len(y_percentile_vals)-1]
#print(X)
#print(y)
X = np.asarray(X).reshape(-1, 1)
y = np.asarray(y).reshape(-1, 1)
degree=80
polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
polyreg.fit(X,y)
plt.figure()
plt.scatter(D3_g_median_mags, D3_sum_intrinsic_rms, s = 0.01)
plt.xlabel('Median <g> Magnitudes')
plt.ylabel('Weighted Summed Intrinsic RMS Squared')
plt.xlim(17, 24.5)
plt.ylim(-0.5, 1)
plt.plot(X,polyreg.predict(X), color="black")
plt.axhline(y = 0.0025, color = "crimson")
plt.axhline(y = 0.01, color = "blue")
plt.scatter(D3_x_percentile_vals, D3_y_percentile_vals, s = 5, c = 'red')
plt.title("Polynomial regression with degree "+str(degree))
plt.show()

# |%%--%%| <tmHlxZxdRA|ZZa9hRNWhQ>

#FOR D4
X = D4_x_percentile_vals#[:len(x_percentile_vals) -1]
y = D4_y_percentile_vals#[:len(y_percentile_vals)-1]
#print(X)
#print(y)
X = np.asarray(X).reshape(-1, 1)
y = np.asarray(y).reshape(-1, 1)
degree=47
polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
polyreg.fit(X,y)
plt.figure()
plt.scatter(D4_g_median_mags, D4_sum_intrinsic_rms, s = 0.01)
plt.xlabel('Median <g> Magnitudes')
plt.ylabel('Weighted Summed Intrinsic RMS Squared')
plt.xlim(17, 24.5)
plt.ylim(-0.5, 1)
plt.plot(X,polyreg.predict(X), color="black")
plt.axhline(y = 0.0025, color = "crimson")
plt.axhline(y = 0.01, color = "blue")
plt.scatter(D4_x_percentile_vals, D4_y_percentile_vals, s = 5, c = 'red')
plt.title("Polynomial regression with degree "+str(degree))
plt.show()

# |%%--%%| <ZZa9hRNWhQ|s2zUkzFDym>

#to find best fit
RMSE_vals = []
degree_vals = []
candidate_nums = []
X = np.asarray(X).reshape(-1, 1)
y = np.asarray(y).reshape(-1, 1)
for degree in range(1,100):
    polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
    polyreg.fit(X,y)
    RMSE_vals.append(sklearn.metrics.mean_squared_error(y,polyreg.predict(X)))
    degree_vals.append(degree)
    #candidate parts
    variable_candidates = []
    for i in range(len(D2_sum_intrinsic_rms)):
        mags_arr = D2_g_median_mags[i]
        mags_arr = mags_arr.reshape(1, -1)
        if(D2_sum_intrinsic_rms[i] > polyreg.predict(mags_arr)):
            variable_candidates.append(i)
    candidate_nums.append(len(variable_candidates))

# |%%--%%| <s2zUkzFDym|j03CTo1eWc>

plt.figure(figsize = (9, 12))
plt.scatter(degree_vals, candidate_nums)
plt.show()

# |%%--%%| <j03CTo1eWc|76PIpYEwRC>
"""°°°
Categorizing variables into subtle, intermediate and highest for each field.
°°°"""
# |%%--%%| <76PIpYEwRC|qstHZE2Dpg>

D1_marginal_variables = []
D1_intermediate_variables = []
D1_extreme_variables = []
D1_cands_over_half = []
for i in range(len(D1_sum_intrinsic_rms)):
    mags_arr = D1_g_median_mags[i]
    mags_arr = mags_arr.reshape(1, -1)
    if(D1_sum_intrinsic_rms[i] > polyreg.predict(mags_arr) and D1_sum_intrinsic_rms[i] <= 0.0025 and D1_g_rms[i] > 0.1 and D1_r_rms[i] > 0.1):
        D1_marginal_variables.append(fields_1[i])
    if(D1_sum_intrinsic_rms[i] > polyreg.predict(mags_arr) and D1_sum_intrinsic_rms[i] > 0.0025 and D1_sum_intrinsic_rms[i] <= 0.01 and D1_g_rms[i]> 0.1 and D1_r_rms[i]>0.1):
        D1_intermediate_variables.append(fields_1[i])
    if(D1_sum_intrinsic_rms[i] > polyreg.predict(mags_arr) and D1_sum_intrinsic_rms[i] > 0.01 and D1_g_rms[i]>0.1 and D1_r_rms[i]>0.1):
        D1_extreme_variables.append(fields_1[i])
    if(D1_sum_intrinsic_rms[i] > 0.5):
        D1_cands_over_half.append(fields_1[i])

# |%%--%%| <qstHZE2Dpg|hAUVg2TPB2>

D2_marginal_variables = []
D2_intermediate_variables = []
D2_extreme_variables = []
D2_cands_over_half = []
for i in range(len(D2_sum_intrinsic_rms)):
    mags_arr = D2_g_median_mags[i]
    mags_arr = mags_arr.reshape(1, -1)
    if(D2_sum_intrinsic_rms[i] > polyreg.predict(mags_arr) and D2_sum_intrinsic_rms[i] <= 0.0025 and D2_g_rms[i]>0.1 and D2_r_rms[i]>0.1):
        D2_marginal_variables.append(fields_2[i])
    if(D2_sum_intrinsic_rms[i] > polyreg.predict(mags_arr) and D2_sum_intrinsic_rms[i] > 0.0025 and D2_sum_intrinsic_rms[i] <= 0.01 and D2_g_rms[i]>0.1 and D2_r_rms[i]>0.1):
        D2_intermediate_variables.append(fields_2[i])
    if(D2_sum_intrinsic_rms[i] > polyreg.predict(mags_arr) and D2_sum_intrinsic_rms[i] > 0.01 and D2_g_rms[i]>0.1 and D2_r_rms[i]>0.1):
        D2_extreme_variables.append(fields_2[i])
    if(D2_sum_intrinsic_rms[i] > 0.5):
        D2_cands_over_half.append(fields_2[i])

# |%%--%%| <hAUVg2TPB2|v3bfmLr9Ia>

D3_marginal_variables = []
D3_intermediate_variables = []
D3_extreme_variables = []
D3_cands_over_half = []
for i in range(len(D3_sum_intrinsic_rms)):
    mags_arr = D3_g_median_mags[i]
    mags_arr = mags_arr.reshape(1, -1)
    if(D3_sum_intrinsic_rms[i] > polyreg.predict(mags_arr) and D3_sum_intrinsic_rms[i] <= 0.0025 and D3_g_rms[i]>0.1 and D3_r_rms[i]>0.1):
        D3_marginal_variables.append(fields_3[i])
    if(D3_sum_intrinsic_rms[i] > polyreg.predict(mags_arr) and D3_sum_intrinsic_rms[i] > 0.0025 and D3_sum_intrinsic_rms[i] <= 0.01 and D3_g_rms[i]>0.1 and D3_r_rms[i]>0.1):
        D3_intermediate_variables.append(fields_3[i])
    if(D3_sum_intrinsic_rms[i] > polyreg.predict(mags_arr) and D3_sum_intrinsic_rms[i] > 0.01 and D3_g_rms[i]>0.1 and D3_r_rms[i]>0.1):
        D3_extreme_variables.append(fields_3[i])
    if(D3_sum_intrinsic_rms[i] > 0.5):
        D3_cands_over_half.append(fields_3[i])

# |%%--%%| <v3bfmLr9Ia|VNkhGnZaEC>

D4_marginal_variables = []
D4_intermediate_variables = []
D4_extreme_variables = []
D4_cands_over_half = []
for i in range(len(D4_sum_intrinsic_rms)):
    mags_arr = D4_g_median_mags[i]
    mags_arr = mags_arr.reshape(1, -1)
    if(D4_sum_intrinsic_rms[i] > polyreg.predict(mags_arr) and D4_sum_intrinsic_rms[i] <= 0.0025 and D4_g_rms[i]>0.1 and D4_r_rms[i]>0.1):
        D4_marginal_variables.append(fields_4[i])
    if(D4_sum_intrinsic_rms[i] > polyreg.predict(mags_arr) and D4_sum_intrinsic_rms[i] > 0.0025 and D4_sum_intrinsic_rms[i] <= 0.01 and D4_g_rms[i]>0.1 and D4_r_rms[i]>0.1):
        D4_intermediate_variables.append(fields_4[i])
    if(D4_sum_intrinsic_rms[i] > polyreg.predict(mags_arr) and D4_sum_intrinsic_rms[i] > 0.01 and D4_g_rms[i]>0.1 and D4_r_rms[i]>0.1):
        D4_extreme_variables.append(fields_4[i])
    if(D4_sum_intrinsic_rms[i] > 0.5):
        D4_cands_over_half.append(fields_4[i])

# |%%--%%| <VNkhGnZaEC|lhIIeFzZo9>

#usually don't download more than 600 at a time or the computer will crash
#reminder: CHARGE THE COMPUTER WHILE DOING THIS
def pdf_unfolded_light_curve(objects, pdf_name):
    """
    Function to create pdfs of unfolded light curve plots for all objects of interest
    
    Parameters
    ---
    objects: array of indicies of objects of interest in file_list
    pdf_name: name of pdf
    
    Returns
    ---
    pdf of magnitude vs MJD for all objects
    """

    plt.ioff()
    with PdfPages(pdf_name) as pdf:
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
            pdf.savefig()
            plt.close()

# |%%--%%| <lhIIeFzZo9|btepmSfYQb>

#Once all the code below has run, USE THIS TO CHECK THE LENGTHS OF EACH SO YOU DON'T CRASH THE COMPUTER WHEN RUNNING THE CODE BLOCK BELOW
print(len(D1_marginal_variables))
print(len(D2_marginal_variables))
print(len(D3_marginal_variables))
print(len(D4_marginal_variables))

print(len(D1_intermediate_variables))
print(len(D2_intermediate_variables))
print(len(D3_intermediate_variables))
print(len(D4_intermediate_variables))

print(len(D1_extreme_variables))
print(len(D2_extreme_variables))
print(len(D3_extreme_variables))
print(len(D4_extreme_variables))

# |%%--%%| <btepmSfYQb|BHh6irkp1s>

#RUN THIS LAST
pdf_unfolded_light_curve([503], "new503.pdf")

# |%%--%%| <BHh6irkp1s|IVKFIdFX6D>

fields[18994]

# |%%--%%| <IVKFIdFX6D|TEXSh55kTg>

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

# |%%--%%| <TEXSh55kTg|B3iBRvi6n7>
"""°°°
Downloading indicies of candidate variables into excel
°°°"""
# |%%--%%| <B3iBRvi6n7|iYZZfb8FYI>

import pandas as pd

# |%%--%%| <iYZZfb8FYI|hNi3zwHn1y>

D1_marginal_df = pd.DataFrame(D1_marginal_variables, columns=['ID'])
D1_intermediate_df = pd.DataFrame(D1_intermediate_variables, columns=['ID'])
D1_extreme_df = pd.DataFrame(D1_extreme_variables, columns=['ID'])

D2_marginal_df = pd.DataFrame(D2_marginal_variables, columns=['ID'])
D2_intermediate_df = pd.DataFrame(D2_intermediate_variables, columns=['ID'])
D2_extreme_df = pd.DataFrame(D2_extreme_variables, columns=['ID'])

D3_marginal_df = pd.DataFrame(D3_marginal_variables, columns=['ID'])
D3_intermediate_df = pd.DataFrame(D3_intermediate_variables, columns=['ID'])
D3_extreme_df = pd.DataFrame(D3_extreme_variables, columns=['ID'])

D4_marginal_df = pd.DataFrame(D4_marginal_variables, columns=['ID'])
D4_intermediate_df = pd.DataFrame(D4_intermediate_variables, columns=['ID'])
D4_extreme_df = pd.DataFrame(D4_extreme_variables, columns=['ID'])

# |%%--%%| <hNi3zwHn1y|Kb2wrWAa8F>

'''
with pd.ExcelWriter('/Users/joannezhao/sip2021/AST-07_UnfoldedLightCurveNotes.xlsx') as writer:
    D1_marginal_df.to_excel(writer, sheet_name = 'D1_marginal')
    D1_intermediate_df.to_excel(writer, sheet_name = 'D1_intermediate')
    D1_extreme_df.to_excel(writer, sheet_name = 'D1_extreme')
    D2_marginal_df.to_excel(writer, sheet_name = 'D2_marginal')
    D2_intermediate_df.to_excel(writer, sheet_name = 'D2_intermediate')
    D2_extreme_df.to_excel(writer, sheet_name = 'D2_extreme')
    D3_marginal_df.to_excel(writer, sheet_name = 'D3_marginal')
    D3_intermediate_df.to_excel(writer, sheet_name = 'D3_intermediate')
    D3_extreme_df.to_excel(writer, sheet_name = 'D3_extreme')
    D4_marginal_df.to_excel(writer, sheet_name = 'D4_marginal')
    D4_intermediate_df.to_excel(writer, sheet_name = 'D4_intermediate')
    D4_extreme_df.to_excel(writer, sheet_name = 'D4_extreme')
'''

# |%%--%%| <Kb2wrWAa8F|WyJN3tGZSO>
"""°°°
The following code creates color color diagrams and histograms.
°°°"""
# |%%--%%| <WyJN3tGZSO|QFCrCxbPAM>

plt.ion()
subtle_variables = []
intermediate_variables = []
highest_variables = []
#initialize:
subtle_variables.extend(D1_subtle_variables)
subtle_variables.extend(D2_subtle_variables)
subtle_variables.extend(D3_subtle_variables)
subtle_variables.extend(D4_subtle_variables)
#initialize:
intermediate_variables.extend(D1_intermediate_variables)
intermediate_variables.extend(D2_intermediate_variables)
intermediate_variables.extend(D3_intermediate_variables)
intermediate_variables.extend(D4_intermediate_variables)
#intitialize:
highest_variables.extend(D1_subtle_variables)
highest_variables.extend(D2_subtle_variables)
highest_variables.extend(D3_subtle_variables)
highest_variables.extend(D4_subtle_variables)
def getx(object_index,string_val):
    obj = all_mags[object_index]
    x = None
    y = None
    if string_val == 'u-g':
        x = np.median(obj[0])
        y = np.median(obj[1])
    elif string_val == 'r-i':
        x = np.median(obj[2])
        y = np.median(obj[3])
    elif string_val == 'i-z':
        x = np.median(obj[3])
        y = np.median(obj[4])
    elif string_val == 'g-r':
        x = np.median(obj[1])
        y = np.median(obj[2])
    val = x - y
    return val
x1_vals = [[],[],[],[]]
y1_vals = [[],[],[],[]]
x2_vals = [[],[],[],[]]
y2_vals = [[],[],[],[]]
x3_vals = [[],[],[],[]]
y3_vals = [[],[],[],[]]
for star_index in range(len(all_mags)):
    if star_index in subtle_variables:
        y1_vals[1].append(getx(star_index,'u-g'))
        x1_vals[1].append(getx(star_index,'g-r'))
        y2_vals[1].append(getx(star_index,'r-i'))
        x2_vals[1].append(getx(star_index,'g-r'))
        y3_vals[1].append(getx(star_index,'i-z'))
        x3_vals[1].append(getx(star_index,'g-r'))
    if star_index in intermediate_variables:
        y1_vals[2].append(getx(star_index,'u-g'))
        x1_vals[2].append(getx(star_index,'g-r'))
        y2_vals[2].append(getx(star_index,'r-i'))
        x2_vals[2].append(getx(star_index,'g-r'))
        y3_vals[2].append(getx(star_index,'i-z'))
        x3_vals[2].append(getx(star_index,'g-r'))
    if star_index in highest_variables:
        y1_vals[3].append(getx(star_index,'u-g'))
        x1_vals[3].append(getx(star_index,'g-r'))
        y2_vals[3].append(getx(star_index,'r-i'))
        x2_vals[3].append(getx(star_index,'g-r'))
        y3_vals[3].append(getx(star_index,'i-z'))
        x3_vals[3].append(getx(star_index,'g-r'))
    else:
        y1_vals[0].append(getx(star_index,'u-g'))
        x1_vals[0].append(getx(star_index,'g-r'))
        y2_vals[0].append(getx(star_index,'r-i'))
        x2_vals[0].append(getx(star_index,'g-r'))
        y3_vals[0].append(getx(star_index,'i-z'))
        x3_vals[0].append(getx(star_index,'g-r'))
fig = plt.figure(figsize = (10, 20))
ax1 = fig.add_subplot(4, 3, 1)
ax2 = fig.add_subplot(4, 3, 2)
ax3 = fig.add_subplot(4, 3, 3)
ax4 = fig.add_subplot(4, 3, 4)
ax5 = fig.add_subplot(4, 3, 5)
ax6 = fig.add_subplot(4, 3, 6)
ax7 = fig.add_subplot(4, 3, 7)
ax8 = fig.add_subplot(4, 3, 8)
ax9 = fig.add_subplot(4, 3, 9)
#spacer --> green for highest, red for subtle, blue for intermediate, black for nonvariables
ax1.scatter(x1_vals[0],y1_vals[0],s = 0.01, c='silver', label = 'nonvariables')
ax1.scatter(x1_vals[2],y1_vals[2],s = 0.05, c='blue', label = 'intermediate')
ax1.scatter(x1_vals[3],y1_vals[3],s = 0.1, c='green', label = 'highest')
ax1.scatter(x1_vals[1],y1_vals[1],s = 0.05, c='crimson', label = 'subtle')
ax2.scatter(x2_vals[0],y2_vals[0],s=0.01,c='silver', label = 'nonvariables')
ax2.scatter(x2_vals[2],y2_vals[2],s=0.05,c='blue', label = 'intermediate')
ax2.scatter(x2_vals[3],y2_vals[3],s=0.1,c='green', label = 'highest')
ax2.scatter(x2_vals[1],y2_vals[1],s=0.05,c='crimson', label = 'subtle')
ax3.scatter(x3_vals[0],y3_vals[0],s=0.01,c='silver', label = 'nonvariables')
ax3.scatter(x3_vals[2],y3_vals[2],s=0.05,c='blue', label = 'intermediate')
ax3.scatter(x3_vals[3],y3_vals[3],s=0.1,c='green', label = 'highest')
ax3.scatter(x3_vals[1],y3_vals[1],s=0.05,c='crimson', label = 'subtle')

ax4.scatter(x1_vals[2],y1_vals[2],s = 0.01, c='blue', label = 'intermediate')
ax4.scatter(x1_vals[3],y1_vals[3],s = 0.05, c='green', label = 'highest')
ax4.scatter(x1_vals[1],y1_vals[1],s = 0.01, c='crimson', label = 'subtle')

ax5.scatter(x2_vals[2],y2_vals[2],s=0.01,c='blue', label = 'intermediate')
ax5.scatter(x2_vals[3],y2_vals[3],s=0.05,c='green', label = 'highest')
ax5.scatter(x2_vals[1],y2_vals[1],s=0.01,c='crimson', label = 'subtle')

ax6.scatter(x3_vals[2],y3_vals[2],s=0.01,c='blue', label = 'intermediate')
ax6.scatter(x3_vals[3],y3_vals[3],s=0.05,c='green', label = 'highest')
ax6.scatter(x3_vals[1],y3_vals[1],s=0.01,c='crimson', label = 'subtle')

x_lower = -0.25
x_higher = 1.5
ax1.set_xlim(x_lower, x_higher)
ax1.set_ylim(-0.5, 3)
ax2.set_xlim(x_lower, x_higher)
ax2.set_ylim(-0.5, 2.5)
ax3.set_xlim(x_lower, x_higher)
ax3.set_ylim(-0.5, 1)
ax4.set_xlim(x_lower, x_higher)
ax4.set_ylim(-0.5, 3)
ax5.set_xlim(x_lower, x_higher)
ax5.set_ylim(-0.5, 2.5)
ax6.set_xlim(x_lower, x_higher)
ax6.set_ylim(-0.5, 1)
ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()
ax5.legend()
ax6.legend()
#ax1.scatter(x1_vals[0,y1_vals[0],c='black')
plt.show()

# |%%--%%| <QFCrCxbPAM|9cDuuTwWhh>


all_cands_over_half = []
all_cands_over_half.extend(D1_cands_over_half)
all_cands_over_half.extend(D1_cands_over_half)
all_cands_over_half.extend(D1_cands_over_half)
all_cands_over_half.extend(D1_cands_over_half)
plot_2_x1 = []
plot_2_y1 = []
plot_2_x2 = []
plot_2_y2 = []
plot_2_x3 = []
plot_2_y3 = []
def getx(object_index,string_val):
    obj = all_mags[object_index]
    x = None
    y = None
    if string_val == 'u-g':
        x = np.median(obj[0])
        y = np.median(obj[1])
    elif string_val == 'r-i':
        x = np.median(obj[2])
        y = np.median(obj[3])
    elif string_val == 'i-z':
        x = np.median(obj[3])
        y = np.median(obj[4])
    elif string_val == 'g-r':
        x = np.median(obj[1])
        y = np.median(obj[2])
    val = x - y
    return val
for i in all_cands_over_half:
    plot_2_y1.append(getx(i,'u-g'))
    plot_2_x1.append(getx(i,'g-r'))
    plot_2_y2.append(getx(i,'r-i'))
    plot_2_x2.append(getx(i,'g-r'))
    plot_2_y3.append(getx(i,'i-z'))
    plot_2_x3.append(getx(i,'g-r'))
ax7.scatter(plot_2_x1,plot_2_y1,s=0.01, c = "black")
ax8.scatter(plot_2_x2,plot_2_y2,s=0.01, c = "black")
ax9.scatter(plot_2_x3,plot_2_y3,s=0.01, c = "black")
ax7.set_xlim(-0.25, 1.5)
ax8.set_xlim(-0.25, 1.5)
ax9.set_xlim(-0.25, 1.5)
ax7.set_ylim(-0.5, 3)
ax8.set_ylim(-0.5, 2.5)
ax9.set_ylim(-0.5, 1)

# |%%--%%| <9cDuuTwWhh|19KGR9lBnW>

not_in = []
for i in f200:
    if i not in variable_candidates:
        print("Not in variable candidates: " + str(i))
        not_in.append(i)

# |%%--%%| <19KGR9lBnW|PVTf3V22CK>

color_color_1 = [i for i in range(len(g_median_mags)) if g_median_mags[i] <= 22.5]
color_color_2 = [i for i in range(len(g_median_mags)) if g_median_mags[i] > 22.5]

# |%%--%%| <PVTf3V22CK|sOMLHvxhFn>

print(len(color_color_1))

# |%%--%%| <sOMLHvxhFn|irLlue0ac3>

print(len(color_color_1))
print(len(color_color_2))

# |%%--%%| <irLlue0ac3|g0hxcLD1VC>

u_median_mags_1 = [np.median(all_mags_2[i][0]) for i in range(len(all_mags_2)) if i in variable_candidates and i in color_color_1]
g_median_mags_1 = [np.median(all_mags_2[i][1]) for i in range(len(all_mags_2)) if i in variable_candidates and i in color_color_1]
r_median_mags_1 = [np.median(all_mags_2[i][2]) for i in range(len(all_mags_2)) if i in variable_candidates and i in color_color_1]
i_median_mags_1 = [np.median(all_mags_2[i][3]) for i in range(len(all_mags_2)) if i in variable_candidates and i in color_color_1]
z_median_mags_1 = [np.median(all_mags_2[i][4]) for i in range(len(all_mags_2)) if i in variable_candidates and i in color_color_1]
u_median_mags_2 = [np.median(all_mags_2[i][0]) for i in range(len(all_mags_2)) if i in variable_candidates and i in color_color_2]
g_median_mags_2 = [np.median(all_mags_2[i][1]) for i in range(len(all_mags_2)) if i in variable_candidates and i in color_color_2]
r_median_mags_2 = [np.median(all_mags_2[i][2]) for i in range(len(all_mags_2)) if i in variable_candidates and i in color_color_2]
i_median_mags_2 = [np.median(all_mags_2[i][3]) for i in range(len(all_mags_2)) if i in variable_candidates and i in color_color_2]
z_median_mags_2 = [np.median(all_mags_2[i][4]) for i in range(len(all_mags_2)) if i in variable_candidates and i in color_color_2]
all_u_median_mags_1 = [np.median(all_mags_2[i][0]) for i in range(len(all_mags_2)) if i in color_color_1]
all_g_median_mags_1 = [np.median(all_mags_2[i][1]) for i in range(len(all_mags_2)) if i in color_color_1]
all_r_median_mags_1 = [np.median(all_mags_2[i][2]) for i in range(len(all_mags_2)) if i in color_color_1]
all_i_median_mags_1 = [np.median(all_mags_2[i][3]) for i in range(len(all_mags_2)) if i in color_color_1]
all_z_median_mags_1 = [np.median(all_mags_2[i][4]) for i in range(len(all_mags_2)) if i in color_color_1]
all_u_median_mags_2 = [np.median(all_mags_2[i][0]) for i in range(len(all_mags_2)) if i in color_color_2]
all_g_median_mags_2 = [np.median(all_mags_2[i][1]) for i in range(len(all_mags_2)) if i in color_color_2]
all_r_median_mags_2 = [np.median(all_mags_2[i][2]) for i in range(len(all_mags_2)) if i in color_color_2]
all_i_median_mags_2 = [np.median(all_mags_2[i][3]) for i in range(len(all_mags_2)) if i in color_color_2]
all_z_median_mags_2 = [np.median(all_mags_2[i][4]) for i in range(len(all_mags_2)) if i in color_color_2]
x_vals_1 = [g_median_mags_1[i] - r_median_mags_1[i] for i in range(len(g_median_mags_1))]
x_vals_2 = [g_median_mags_2[i] - r_median_mags_2[i] for i in range(len(g_median_mags_2))]
u_minus_g_1 = [u_median_mags_1[i] - g_median_mags_1[i] for i in range(len(g_median_mags_1))]
u_minus_g_2 = [u_median_mags_2[i] - g_median_mags_2[i] for i in range(len(g_median_mags_2))]
r_minus_i_1 = [r_median_mags_1[i] - i_median_mags_1[i] for i in range(len(r_median_mags_1))]
r_minus_i_2 = [r_median_mags_2[i] - i_median_mags_2[i] for i in range(len(r_median_mags_2))]
i_minus_z_1 = [i_median_mags_1[i] - z_median_mags_1[i] for i in range(len(i_median_mags_1))]
i_minus_z_2 = [i_median_mags_2[i] - z_median_mags_2[i] for i in range(len(i_median_mags_2))]
x_vals_3 = [all_g_median_mags_1[i] - all_r_median_mags_1[i] for i in range(len(all_g_median_mags_1))]
x_vals_4 = [all_g_median_mags_2[i] - all_r_median_mags_2[i] for i in range(len(all_g_median_mags_2))]
all_u_minus_g_1 = [all_u_median_mags_1[i] - all_g_median_mags_1[i] for i in range(len(all_g_median_mags_1))]
all_u_minus_g_2 = [all_u_median_mags_2[i] - all_g_median_mags_2[i] for i in range(len(all_g_median_mags_2))]
all_r_minus_i_1 = [all_r_median_mags_1[i] - all_i_median_mags_1[i] for i in range(len(all_r_median_mags_1))]
all_r_minus_i_2 = [all_r_median_mags_2[i] - all_i_median_mags_2[i] for i in range(len(all_r_median_mags_2))]
all_i_minus_z_1 = [all_i_median_mags_1[i] - all_z_median_mags_1[i] for i in range(len(all_i_median_mags_1))]
all_i_minus_z_2 = [all_i_median_mags_2[i] - all_z_median_mags_2[i] for i in range(len(all_i_median_mags_2))]

# |%%--%%| <g0hxcLD1VC|WnsWXvoCzN>

fig = plt.figure(figsize = (10, 20))
ax1 = fig.add_subplot(4, 3, 1)
ax2 = fig.add_subplot(4, 3, 2)
ax3 = fig.add_subplot(4, 3, 3)
ax4 = fig.add_subplot(4, 3, 4)
ax5 = fig.add_subplot(4, 3, 5)
ax6 = fig.add_subplot(4, 3, 6)
ax7 = fig.add_subplot(4, 3, 7)
ax8 = fig.add_subplot(4, 3, 8)
ax9 = fig.add_subplot(4, 3, 9)
ax10 = fig.add_subplot(4, 3, 10)
ax11 = fig.add_subplot(4, 3, 11)
ax12 = fig.add_subplot(4, 3, 12)
ax1.scatter(x_vals_1,u_minus_g_1,s = 0.01)
ax2.scatter(x_vals_1,r_minus_i_1,s = 0.01)
ax3.scatter(x_vals_1,i_minus_z_1,s = 0.01)
ax4.scatter(x_vals_2, u_minus_g_2, s = 0.01)
ax5.scatter(x_vals_2, r_minus_i_2, s = 0.01)
ax6.scatter(x_vals_2, i_minus_z_2, s = 0.01)
ax7.scatter(x_vals_3,all_u_minus_g_1,s = 0.01)
ax8.scatter(x_vals_3,all_r_minus_i_1,s = 0.01)
ax9.scatter(x_vals_3,all_i_minus_z_1,s = 0.01)
ax10.scatter(x_vals_4, all_u_minus_g_2, s = 0.01)
ax11.scatter(x_vals_4, all_r_minus_i_2, s = 0.01)
ax12.scatter(x_vals_4, all_i_minus_z_2, s = 0.01)
ax1.set_xlim(-0.5,1.5)
ax4.set_xlim(-0.5, 1.5)
ax7.set_xlim(-0.5,1.5)
ax10.set_xlim(-0.5, 1.5)
ax2.set_xlim(-0.5,1.5)
ax5.set_xlim(-0.5, 1.5)
ax8.set_xlim(-0.5,1.5)
ax11.set_xlim(-0.5, 1.5)
ax3.set_xlim(-0.5,1.5)
ax6.set_xlim(-0.5, 1.5)
ax9.set_xlim(-0.5,1.5)
ax12.set_xlim(-0.5, 1.5)
ax1.set_ylim(-0.3,2.5)
ax4.set_ylim(-0.3,2.5)
ax7.set_ylim(-0.3,2.5)
ax10.set_ylim(-0.3,2.5)
ax2.set_xlim(-0.5,2)
ax5.set_xlim(-0.5,2)
ax8.set_xlim(-0.5,2)
ax11.set_xlim(-0.5,2)
ax3.set_ylim(-0.4,0.9)
ax6.set_ylim(-0.4,0.9)
ax9.set_ylim(-0.4,0.9)
ax12.set_ylim(-0.4,0.9)
ax1.set_xlabel('<g> - <r>')
ax2.set_xlabel('<g> - <r>')
ax3.set_xlabel('<g> - <r>')
ax4.set_xlabel('<g> - <r>')
ax5.set_xlabel('<g> - <r>')
ax6.set_xlabel('<g> - <r>')
ax1.set_ylabel('<u> - <g>')
ax2.set_ylabel('<r> - <i>')
ax3.set_ylabel('<i> - <z>')
ax4.set_ylabel('<u> - <g>')
ax5.set_ylabel('<r> - <i>')
ax6.set_ylabel('<i> - <z>')
ax7.set_xlabel('<g> - <r>')
ax8.set_xlabel('<g> - <r>')
ax9.set_xlabel('<g> - <r>')
ax10.set_xlabel('<g> - <r>')
ax11.set_xlabel('<g> - <r>')
ax12.set_xlabel('<g> - <r>')
ax7.set_ylabel('<u> - <g>')
ax8.set_ylabel('<r> - <i>')
ax9.set_ylabel('<i> - <z>')
ax10.set_ylabel('<u> - <g>')
ax11.set_ylabel('<r> - <i>')
ax12.set_ylabel('<i> - <z>')
plt.show()

# |%%--%%| <WnsWXvoCzN|kv7GDoqkDh>

all_bands_intrinsic_rms_squared = [u_intrinsic_rms_squared,g_intrinsic_rms_squared,r_intrinsic_rms_squared,i_intrinsic_rms_squared,z_intrinsic_rms_squared]
range_vals = np.arange(min(g_median_mags),max(g_median_mags),0.5)
in_bin_medians_g = []
in_bin_intrinsic = []
for i in range(len(range_vals)-1):
    temp_medians = []
    temp_intrinsic = []
    for j in range(len(g_median_mags)):
        if (range_vals[i] < g_median_mags[j] < range_vals[i+1]):
            temp_medians.append(g_median_mags[j])
            temp_intrinsic.append(sum([all_bands_intrinsic_rms_squared[z][j] for z in range(5)]))
    in_bin_medians_g.append(temp_medians)
    in_bin_intrinsic.append(temp_intrinsic)

# |%%--%%| <kv7GDoqkDh|rW8btXvSiR>

def make_histogram(intrinsic_rms_squared,i):
    x_vals = [intrinsic_rms_squared[i][0][0] for i in range(len(intrinsic_rms_squared))]
    plt.figure()
    print('max value: ' + str(max(x_vals)))
    plt.title("Bin: " + str(i))
    plt.xlabel("Summed Intrinsic Rms\u00b2")
    plt.ylabel("Frequency")
    plt.hist(x_vals,bins=300)
    plt.show()

# |%%--%%| <rW8btXvSiR|6E5wYh15fx>

for i in range(len(in_bin_intrinsic)):
    make_histogram(in_bin_intrinsic[i],i)

# |%%--%%| <6E5wYh15fx|CziGoSmmBg>

def plot_rms_mags (band, mags, mags_for_median):
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
    plt.scatter(medians,np.log10(all_rms), s = 0.05)
    plt.scatter(in_bin_medians,np.log10(red_rms), s = 13, c = "crimson")#c = for_color, cmap = 'PuRd'
    if band == 0:
        plt.ylabel(r'$\sigma_%s$' % 'u')
        plt.xlabel(r'$<%s>$' % 'u')
        if mags == D1_modified:
            X = (in_bin_medians[6:len(in_bin_medians)-3])
            y = np.log10(red_rms[6:len(red_rms) - 3])
        elif mags == D2_modified:
            X = (in_bin_medians[6:len(in_bin_medians)-5])
            y = np.log10(red_rms[6:len(red_rms) - 5])
        elif mags == D3_modified:
            X = (in_bin_medians[10:len(in_bin_medians)-6])
            y = np.log10(red_rms[10:len(red_rms) - 6])
        elif mags == D4_modified:
            X = (in_bin_medians[8:len(in_bin_medians)-6])
            y = np.log10(red_rms[8:len(red_rms) - 6])
    elif band == 1:
        plt.ylabel(r'$\sigma_%s$' % 'g')
        plt.xlabel(r'$<%s>$' % 'g')
        if mags == D1_modified:
            X = (in_bin_medians[4:len(in_bin_medians)-2])
            y = np.log10(red_rms[4:len(red_rms) - 2])
        elif mags == D2_modified:
            X = (in_bin_medians[8:len(in_bin_medians)-2])
            y = np.log10(red_rms[8:len(red_rms) - 2])
        elif mags == D3_modified:
            X = (in_bin_medians[9:len(in_bin_medians)-1])
            y = np.log10(red_rms[9:len(red_rms)-1])
        elif mags == D4_modified:
            X = (in_bin_medians[0:len(in_bin_medians)-8])
            y = np.log10(red_rms[0:len(red_rms)-8])
    elif band == 2:
        plt.ylabel(r'$\sigma_%s$' % 'r')
        plt.xlabel(r'$<%s>$' % 'r')
        if mags == D1_modified:
            X = (in_bin_medians[9:len(in_bin_medians)-4])
            y = np.log10(red_rms[9:len(red_rms)-4])
        elif mags == D2_modified:
            X = (in_bin_medians[9:len(in_bin_medians)-2])
            y = np.log10(red_rms[9:len(red_rms)-2])
        elif mags == D3_modified:
            X = (in_bin_medians[12:len(in_bin_medians)-2])
            y = np.log10(red_rms[12:len(red_rms)-2])
        elif mags == D4_modified:
            X = (in_bin_medians[17:len(in_bin_medians)-5])
            y = np.log10(red_rms[17:len(red_rms)-5])
    elif band == 3:
        plt.ylabel(r'$\sigma_%s$' % 'i')
        plt.xlabel(r'$<%s>$' % 'i')
        if mags == D1_modified:
            X = (in_bin_medians[12:len(in_bin_medians)-3])
            y = np.log10(red_rms[12:len(red_rms)-3])
        elif mags == D2_modified:
            X = (in_bin_medians[14:len(in_bin_medians)])
            y = np.log10(red_rms[14:len(red_rms)])
        elif mags == D3_modified:
            X = (in_bin_medians[25:len(in_bin_medians)-1])
            y = np.log10(red_rms[25:len(red_rms)-1])
        elif mags == D4_modified:
            X = (in_bin_medians[16:len(in_bin_medians)-5])
            y = np.log10(red_rms[16:len(red_rms)-5])
    elif band == 4:
        plt.ylabel(r'$\sigma_%s$' % 'z')
        plt.xlabel(r'$<%s>$' % 'z')
        if mags == D1_modified:
            X = (in_bin_medians[10:len(in_bin_medians)-6])
            y = np.log10(red_rms[10:len(red_rms)-6])
        elif mags == D2_modified:
            X = (in_bin_medians[15:len(in_bin_medians)-4])
            y = np.log10(red_rms[15:len(red_rms)-4])
        elif mags == D3_modified:
            X = (in_bin_medians[14:len(in_bin_medians)-4])
            y = np.log10(red_rms[14:len(red_rms)-4])
        elif mags == D4_modified:
            X = (in_bin_medians[14:len(in_bin_medians)-4])
            y = np.log10(red_rms[14:len(red_rms)-4])
    
    X = np.asarray(X).reshape(-1,1)
    y = np.asarray(y).reshape(-1,1)
    degree=4
    polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
    polyreg.fit(X,y)
    plt.plot(X,polyreg.predict(X),color="black")
    plt.title("Polynomial regression with degree "+str(degree))
    #degree=61
    #polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
    #polyreg.fit(X,y)
    #plt.plot(X,polyreg.predict(X),color="black")
    #plt.title("Polynomial regression with degree "+str(degree))
    #print("RMSE: " + str(sklearn.metrics.mean_squared_error(y,polyreg.predict(X))))
    plt.show()
    return medians, all_rms

# |%%--%%| <CziGoSmmBg|alKhn0dmpD>

all_mags

# |%%--%%| <alKhn0dmpD|hWfLXvrHsQ>



# |%%--%%| <hWfLXvrHsQ|bb9ctLyRYp>

def modify_fields(Dx):
    Dx_modified = []
    for i in Dx:
        Dx_modified.append(all_mags[i])
    return Dx_modified
        

# |%%--%%| <bb9ctLyRYp|Iy3a6tb16n>

D1_modified = modify_fields(D1)
D2_modified = modify_fields(D2)
D3_modified = modify_fields(D3)
D4_modified = modify_fields(D4)

# |%%--%%| <Iy3a6tb16n|wKZvRXXjaD>

D1_modified[0][0]

# |%%--%%| <wKZvRXXjaD|6BblRgv070>



# |%%--%%| <6BblRgv070|gHJvSrfybr>

D1_u_medians, D1_u_rms = plot_rms_mags(0, D1_modified, D1_sigma)

# |%%--%%| <gHJvSrfybr|CFylVOdsfh>

def get_error(band,mags,mags_for_median):
    
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
    plt.scatter(medians,np.log10(all_rms), s = 0.05)
    plt.scatter(in_bin_medians,np.log10(red_rms), s = 13, c = "crimson")#c = for_color, cmap = 'PuRd'
    if band == 0:
        plt.ylabel(r'$\sigma_%s$' % 'u')
        plt.xlabel(r'$<%s>$' % 'u')
        if mags == D1_modified:
            X = (in_bin_medians[6:len(in_bin_medians)-3])
            y = np.log10(red_rms[6:len(red_rms) - 3])
        elif mags == D2_modified:
            X = (in_bin_medians[6:len(in_bin_medians)-5])
            y = np.log10(red_rms[6:len(red_rms) - 5])
        elif mags == D3_modified:
            X = (in_bin_medians[10:len(in_bin_medians)-6])
            y = np.log10(red_rms[10:len(red_rms) - 6])
        elif mags == D4_modified:
            X = (in_bin_medians[8:len(in_bin_medians)-6])
            y = np.log10(red_rms[8:len(red_rms) - 6])
    elif band == 1:
        plt.ylabel(r'$\sigma_%s$' % 'g')
        plt.xlabel(r'$<%s>$' % 'g')
        if mags == D1_modified:
            X = (in_bin_medians[4:len(in_bin_medians)-2])
            y = np.log10(red_rms[4:len(red_rms) - 2])
        elif mags == D2_modified:
            X = (in_bin_medians[8:len(in_bin_medians)-2])
            y = np.log10(red_rms[8:len(red_rms) - 2])
        elif mags == D3_modified:
            X = (in_bin_medians[9:len(in_bin_medians)-1])
            y = np.log10(red_rms[9:len(red_rms)-1])
        elif mags == D4_modified:
            X = (in_bin_medians[0:len(in_bin_medians)-8])
            y = np.log10(red_rms[0:len(red_rms)-8])
    elif band == 2:
        plt.ylabel(r'$\sigma_%s$' % 'r')
        plt.xlabel(r'$<%s>$' % 'r')
        if mags == D1_modified:
            X = (in_bin_medians[9:len(in_bin_medians)-4])
            y = np.log10(red_rms[9:len(red_rms)-4])
        elif mags == D2_modified:
            X = (in_bin_medians[9:len(in_bin_medians)-2])
            y = np.log10(red_rms[9:len(red_rms)-2])
        elif mags == D3_modified:
            X = (in_bin_medians[12:len(in_bin_medians)-2])
            y = np.log10(red_rms[12:len(red_rms)-2])
        elif mags == D4_modified:
            X = (in_bin_medians[17:len(in_bin_medians)-5])
            y = np.log10(red_rms[17:len(red_rms)-5])
    elif band == 3:
        plt.ylabel(r'$\sigma_%s$' % 'i')
        plt.xlabel(r'$<%s>$' % 'i')
        if mags == D1_modified:
            X = (in_bin_medians[12:len(in_bin_medians)-3])
            y = np.log10(red_rms[12:len(red_rms)-3])
        elif mags == D2_modified:
            X = (in_bin_medians[14:len(in_bin_medians)])
            y = np.log10(red_rms[14:len(red_rms)])
        elif mags == D3_modified:
            X = (in_bin_medians[25:len(in_bin_medians)-1])
            y = np.log10(red_rms[25:len(red_rms)-1])
        elif mags == D4_modified:
            X = (in_bin_medians[16:len(in_bin_medians)-5])
            y = np.log10(red_rms[16:len(red_rms)-5])
    elif band == 4:
        plt.ylabel(r'$\sigma_%s$' % 'z')
        plt.xlabel(r'$<%s>$' % 'z')
        if mags == D1_modified:
            X = (in_bin_medians[10:len(in_bin_medians)-6])
            y = np.log10(red_rms[10:len(red_rms)-6])
        elif mags == D2_modified:
            X = (in_bin_medians[15:len(in_bin_medians)-4])
            y = np.log10(red_rms[15:len(red_rms)-4])
        elif mags == D3_modified:
            X = (in_bin_medians[14:len(in_bin_medians)-4])
            y = np.log10(red_rms[14:len(red_rms)-4])
        elif mags == D4_modified:
            X = (in_bin_medians[14:len(in_bin_medians)-4])
            y = np.log10(red_rms[14:len(red_rms)-4])
    
    X = np.asarray(X).reshape(-1,1)
    y = np.asarray(y).reshape(-1,1)
    degree=4
    polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
    polyreg.fit(X,y)
  
    plt.plot(X,polyreg.predict(X),color="black")
    plt.title("Polynomial regression with degree "+str(degree))
    



    plt.show()
    #return medians, all_rms
  
    #Optional: get medians
    medians_obj_band = all_mags.copy()
    for i in range(len(all_mags)):
        for j in range(len(all_mags[i])):
            medians_obj_band[i][j] = np.median(all_mags[i][j])
            
    Y_X = None
    error_measurements = medians_obj_band.copy()
    for i in range(len(medians_obj_band)):
        for j in range(len(medians_obj_band[i])):
            Y_X = polyreg.predict(np.asarray(medians_obj_band[i][j]).reshape(-1,1))
            error_measurements[i][j] = Y_X
            
            
    error_measurements_list = all_mags.copy()
    for i in range(len(error_measurements_list)):
        for j in range(5):
            error_measurements_list[i][j] = error_measurements[i][j][0][0]

        
    return error_measurements_list
    

# |%%--%%| <CFylVOdsfh|ezIe9tOLM9>

def get_error_by_field(field_num):
    if field_num == 1:
        u_errors = get_error(0,D1_modified,D1_sigma)
        g_errors = get_error(1,D1_modified,D1_sigma)
        r_errors = get_error(2,D1_modified,D1_sigma)
        i_errors = get_error(3,D1_modified,D1_sigma)
        z_errors = get_error(4,D1_modified,D1_sigma)
    elif field_num == 2:
        u_errors = get_error(0,D2_modified,D2_sigma)
        g_errors = get_error(1,D2_modified,D2_sigma)
        r_errors = get_error(2,D2_modified,D2_sigma)
        i_errors = get_error(3,D2_modified,D2_sigma)
        z_errors = get_error(4,D2_modified,D2_sigma)
        
    elif field_num == 3:
        u_errors = get_error(0,D3_modified,D3_sigma)
        g_errors = get_error(1,D3_modified,D3_sigma)
        r_errors = get_error(2,D3_modified,D3_sigma)
        i_errors = get_error(3,D3_modified,D3_sigma)
        z_errors = get_error(4,D3_modified,D3_sigma)
    elif field_num == 4:
        u_errors = get_error(0,D4_modified,D1_sigma)
        g_errors = get_error(1,D4_modified,D4_sigma)
        r_errors = get_error(2,D4_modified,D4_sigma)
        i_errors = get_error(3,D4_modified,D4_sigma)
        z_errors = get_error(4,D4_modified,D4_sigma)
    else:
        print("Error!")
        
    all_errors = [u_errors,g_errors,r_errors,i_errors,z_errors]
    return all_errors

# |%%--%%| <ezIe9tOLM9|ot9LXgJr7C>

def get_error_by_object(all_errors_1,all_errors_2,all_errors_3,all_errors_4):
    all_errors = all_mags.copy()
    #CHECK: all_errors.clear()
    counter_1 = 0
    counter_2 = 0
    counter_3 = 0
    counter_4 = 0
    for i in range(len(all_mags)):
        if i in D1:
            for j in range(len(5)):
                all_errors[i][j] = all_errors_1[j][counter_1]
            counter_1 += 1
        elif i in D2:
            for j in range(len(5)):
                all_errors[i][j] = all_errors_2[j][counter_2]
            counter_2 += 1
        elif i in D3:
            for j in range(len(5)):
                all_errors[i][j] = all_errors_3[j][counter_3]
            counter_3 += 1
        elif i in D4:
            for j in range(len(5)):
                all_errors[i][j] = all_errors_4[j][counter_4]
            counter_4 += 1
        else:
            print("Error!")
            
        return all_errors
        

# |%%--%%| <ot9LXgJr7C|k4ezMIZkQj>

all_errors_field_1 = get_error_by_field(1)
all_errors_field_2 = get_error_by_field(2)
all_errors_field_3 = get_error_by_field(3)
all_errors_field_4 = get_error_by_field(4)


# |%%--%%| <k4ezMIZkQj|nSR7YNYSTN>

get_error_by_object(all_errors_field_1,all_errors_field_2,all_errors_field_3,all_errors_field_4)

# |%%--%%| <nSR7YNYSTN|7wbwAGSFhk>

D1[0:4]

# |%%--%%| <7wbwAGSFhk|ciMI6ca2EJ>

D1_

# |%%--%%| <ciMI6ca2EJ|gbYUcRYnqX>



# |%%--%%| <gbYUcRYnqX|vI8nIFoXP1>

t, mags, dy, filts = data_format(26059)

#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib as mpl
#mpl.style.use('ggplot')
#mpl.rc('axes', color_cycle=["#4C72B0", "#55A868", "#C44E52","#8172B2", "#CCB974"]) colors for each band
#from gatspy.periodic import LombScargleMultiband
#from gatspy import datasets, periodic

# Choose a Sesar 2010 object to base our fits on
#lcid = 1019544
#rrlyrae = datasets.RRLyraeGenerated(lcid, random_state=0)

# Generate data in a 6-month observing season
#Nobs = 60
#rng = np.random.RandomState(0)

#nights = np.arange(180)
#rng.shuffle(nights)
#nights = nights[:Nobs]

#t = 57000 + nights + 0.05 * rng.randn(Nobs)
#dy = 0.06 + 0.01 * rng.randn(Nobs)
#mags = np.array([rrlyrae.generated(band, t, err=dy, corrected=False)for band in 'ugriz'])
"""
All the above are not necessary for the NGVS_legacy work, becuase you do have realistic data.
I generate mock datesets just to test to code and show an example result.
"""

"""
Here start the fitting process. 
t, mags, dy, filts are all N-d arraies of your observation data. They mean observation time in MJD, apperent magnitudes, errors, and filter list, respectively
If you do not have error list, set dy to be all ones.
filts could be anything in the format of np.array(['u','u','g', ..., 'z']), as long as its elements are filter names and its length equal to the mags array
"""
#example of appropriate format
#filts = np.take(list('ugriz'), np.arange(Nobs), mode='wrap') 
# 
#mags = mags[np.arange(Nobs) % 5, np.arange(Nobs)]
#masks = [(filts == band) for band in 'ugriz']#separate ugriz to 5 sublists

#below is the necessary code
periods = np.linspace(0.1, 1.0, 100000) # This defines the search range of your period, you can specify it at your will. These are in days.


#2 different ways of fitting
##model = periodic.NaiveMultiband(BaseModel=periodic.LombScargleFast) 
# specify the method to be the naive multiband LS, which means you fit data in each band separately, and get a score list for each band.
# serves as a good first try on your NGVS_legacy data
#above is the fastest way. x axis is period, y axis is score. 5d array output. real variable should have same peak for all bands
##model.fit(t, mags, dy, filts) 
##P = model.scores(periods) 
# This is the fitting score list you want. 
# It is a 5xN array, where N is number of periods tried, P[0] is the fit socres of periods with your u band data. And so on for P[1] for g, P[2] for i, ...
# Each element in P[i] correspond to a period in the array periods you input. The closer to 1, the better.


#all bands done together
LS_multi = periodic.LombScargleMultiband(Nterms_base=1, Nterms_band=0)#initiate structure variable
LS_multi.fit(t, mags, dy, filts)#input our data
P_multi = LS_multi.periodogram(periods)#function where input is periods
# A non-naive way of multiband fitting. This time all data from all bands are fitted simultaneously, means you do not get scores for each band separately.
# P_multi will be a 1-d array, has equal length to your input periods. The maximum value in P_multi will be the best fit, and its corresponding period will be the best period.

#do both!!

"""
From here are visualization of the results.

fig = plt.figure(figsize=(10, 4))
gs = plt.GridSpec(5, 2, left=0.07, right=0.95, bottom=0.15,
                  wspace=0.1, hspace=0.6)
ax = [fig.add_subplot(gs[:, 0]),
      fig.add_subplot(gs[:-2, 1]),
      fig.add_subplot(gs[-2:, 1])]

for band, mask in zip('ugriz', masks):
    ax[0].errorbar((t[mask] / rrlyrae.period) % 1, mags[mask], dy[mask],
                   fmt='.', label=band)
ax[0].set_ylim(18, 14.5)
ax[0].legend(loc='upper left', fontsize=12, ncol=3)
ax[0].set_title('Folded Data, 1 band per night (P={0:.3f} days)'
                ''.format(rrlyrae.period), fontsize=12)
ax[0].set_xlabel('phase')
ax[0].set_ylabel('magnitude')

for i, band in enumerate('ugriz'):
    offset = 4 - i
    ax[1].plot(periods, P[band] + offset, lw=1)
    ax[1].text(0.89, 1 + offset, band, fontsize=10, ha='right', va='top')
ax[1].set_title('Standard Periodogram in Each Band', fontsize=12)
ax[1].yaxis.set_major_formatter(plt.NullFormatter())
ax[1].xaxis.set_major_formatter(plt.NullFormatter())
ax[1].set_ylabel('power + offset')


ax[2].plot(periods, P_multi, lw=1, color='gray')

ax[2].set_title('Multiband Periodogram', fontsize=12)
ax[2].set_yticks([0, 0.5, 1.0])
ax[2].set_ylim(0, 1.0)
ax[2].yaxis.set_major_formatter(plt.NullFormatter())
ax[2].set_xlabel('Period (days)')
ax[2].set_ylabel('power')
"""


plt.ion() #turns figure display on
plt.figure()
plt.scatter(periods, P_multi, s = 0.05)

# |%%--%%| <vI8nIFoXP1|vGxgW8eX08>

best_period = max(P_multi)
for i in range(len(P_multi)):
    if P_multi[i] == best_period:
        index = i
print(periods[index])

# |%%--%%| <vGxgW8eX08|v0WYSyk8b7>

def not_pdf_folded_light_curve(obj, period):
    #rr lyrae has variation in g band of 0.6 or 0.5 in unfolded light curve
    plt.ion()
    #with PdfPages(pdf_name) as pdf:
        #for star_object in objects:
    plt.figure(figsize = (9, 12))
    plt.xlabel('Phase')
    plt.ylabel('Magnitude')
    plt.gca().invert_yaxis()
    peak_mag = min(all_mags[obj][0])
    u_phase = [(i%period)/period for i in all_mjd[obj][0]]
    g_phase = [(i%period)/period for i in all_mjd[obj][1]]
    r_phase = [(i%period)/period for i in all_mjd[obj][2]]
    i_phase = [(i%period)/period for i in all_mjd[obj][3]]
    z_phase = [(i%period)/period for i in all_mjd[obj][4]]
    for i in range(len(all_mags[obj][0])):
        if all_mags[obj][0][i] == peak_mag:
            peak_index = i
    u_phase_change = [u_phase[i]- u_phase[peak_index] if u_phase[i]- u_phase[peak_index] >= 0 else u_phase[i]- u_phase[peak_index] + 1 for i in range(len(u_phase))]
    g_phase_change = [g_phase[i]- u_phase[peak_index] if g_phase[i]- u_phase[peak_index] >= 0 else g_phase[i]- u_phase[peak_index] + 1 for i in range(len(g_phase))]
    r_phase_change = [r_phase[i]- u_phase[peak_index] if r_phase[i]- u_phase[peak_index] >= 0 else r_phase[i]- u_phase[peak_index] + 1 for i in range(len(r_phase))]
    i_phase_change = [i_phase[i]- u_phase[peak_index] if i_phase[i]- u_phase[peak_index] >= 0 else i_phase[i]- u_phase[peak_index] + 1 for i in range(len(i_phase))]
    z_phase_change = [z_phase[i]- u_phase[peak_index] if z_phase[i]- u_phase[peak_index] >= 0 else z_phase[i]- u_phase[peak_index] + 1 for i in range(len(z_phase))]
        #f*p^2/length of observations, f is about 0.1, f = phase error for individual cycle
    plt.scatter(u_phase_change, all_mags[obj][0], s = 5, c = 'blue', label = 'u')
    plt.scatter(g_phase_change, all_mags[obj][1], s = 5, c = 'green', label = 'g')
    plt.scatter(r_phase_change, all_mags[obj][2], s = 5, c = 'purple', label = 'r')
    plt.scatter(i_phase_change, all_mags[obj][3], s = 5, c = 'gold', label = 'i')
    plt.scatter(z_phase_change, all_mags[obj][4], s = 5, c = 'tab:red', label = 'z')
    plt.legend()
    plt.title(fnames[obj] + " (" + str(obj) + ")")
        #pdf.savefig()
        #plt.close()
    #0.6 seconds ideal step
    # i/p_true + n where n is an integer all reciprocated is the beat frequency

# |%%--%%| <v0WYSyk8b7|Gqg1dsi14a>

not_pdf_folded_light_curve(26059, periods[index])

# |%%--%%| <Gqg1dsi14a|VF1boHmBji>

def lomb_scargle_list (object_list):
max_over_mean = []
for cur_object in object_list:
t, mags, dy, filts = data_format(cur_object)
periods = np.linspace(0.1, 1.0, 100000)
LS_multi = periodic.LombScargleMultiband(Nterms_base=1, Nterms_band=0)#initiate structure variable
LS_multi.fit(t, mags, dy, filts)#input our data
P_multi = LS_multi.periodogram(periods)#function where input is periods
max_over_mean.append(np.max(P_multi)/np.mean(P_multi))
return max_over_mean 
