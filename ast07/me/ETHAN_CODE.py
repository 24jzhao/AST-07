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








#%matplotlib notebook

full_file_list = os.listdir('../full') # creates a Python list containing the file paths for every object's time series

# |%%--%%| <8X5ZoLSZ5N|YN09D7F4yY>

df = pd.read_csv('all_categ_vals.csv')
#df.drop('Unnamed: 0',inplace=True,axis=1)
all_categ_vals = df.values.tolist()

#Get indicies for value set
three_object_indexes = [w for w,x,y,z in all_categ_vals if (x==3 and y == 3 and z == 3)]


#print(three_object_indexes)
#len(three_object_indexes)

# |%%--%%| <YN09D7F4yY|I219QxkTi5>



#|%%--%%| <I219QxkTi5|WnYqqsnf36>



