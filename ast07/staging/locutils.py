import numpy
from numpy.random import default_rng
import pandas as pd

def randomgen(rval, rsize):
    rng = default_rng()
    rint = rng.integers(low=0, high=len(rval), size=rsize)
    return rint
