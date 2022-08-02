"""
Python script to generate a large dataset for evolution data plotting
(see PNAS paper fig 3 for reference).

CuPy (& CUDA toolkit) required.
"""

import netCDF4 as nc
import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
from utils_cupy import *  # utils_cupy functions are used - see script
                          # docstring for more information

def bulk_compute(years, T):
    """
    Bulk computes all data for every month-offset (0 to 11) of every year.

    Parameters:
    years:      NumPy array of all the years for which data will be generated
    T:          a 3-dimensional time-series data CuPy array containing the 
                daily temperature data for every grid point. Realistically
                this would need to be T from temp_data_1948_2021.npy
                as per generated from nino.py. This is because year_series
                will be called for all the years computed and the source
                T array will need to cover all the years since 1948 to stay 
                consistent with the parameter settings in year_series (see
                utils.py).
    
    Returns:
    results:    NumPy array of the computation results. Dimensions:
                (no. years, no. of monthly offsets, no. of result types,
                no. of grids outside the el-nino basin). 
                The result types in this case are the in-weights and
                number of in-links for each grid (or node) outside the el-
                nino basin.
    """

    # test fetching results from year_series w/ year 1972
    print("Initialising dimensions and checking for errors...")
    C, T_in, T_out = year_series(T, lat, lon, 1972)  

    # monthly offsets by the number of days
    month_days = [0, 31, 28, 31, 30, 31, 
                  30, 31, 31, 30, 31, 30]  # no leap years (consistent)

    # initialise an empty results array (2077 is arbitrary - just for
    # ease of indexing years when plotting)
    results = np.empty([2077, 
                        12, 2, cp.shape(T_out)[0]])
    print("The empty result array has shape: ",
           cp.shape(results))

    # computing data
    for shift_index in range(12):  # monthly offset from 0 to 11
        """
        instead of messing with the time indicies we simply roll 
        the data T backwards by an amount of days, so the same 
        time indicies year_series is working off of now refers to
        the dates of a given year offset by some amount of months.
        e.g. a year in the original data would go from June to June,
        but a 2-month offset makes it August to August, and so on.
        """
        T = cp.roll(T, -month_days[shift_index], axis = 0)  # bit slow

        for year in years:
            print("Computing year:", year, "; shift:", shift_index)

            # get results for a year
            C, T_in, T_out = year_series(T, lat, lon, year) 

            # calculate in-weights and in-links
            in_N = cp.sum((C[:,:,1] < 151), axis = 0)
            in_C = cp.sum(C[:,:,0], axis = 0)

            # store results in corresponding year and offset amount
            results[year, shift_index, 0, :] = cp.asnumpy(in_C)
            results[year, shift_index, 1, :] = cp.asnumpy(in_N)
    
    return results  # note that we are returning NumPy not CuPy

# loading the data generate in nino.py
with open('temp_data_1948_2021.npy', 'rb') as f:
    T = np.load(f)
    lat = np.load(f)
    lon = np.load(f)

# convert to CuPy
T = cp.asarray(T)
lat = cp.asarray(lat)
lon = cp.asarray(lon)

# we only need this interval to plot the evolution graph
years_to_compute = np.arange(1949, 2017)
results = bulk_compute(years_to_compute, T)
# these might not be relevant but I've kept them for now
lat = cp.asnumpy(lat)
lon = cp.asnumpy(lon)

# save in npy file for plotting and more
with open('bulk_temp_1949_2016.npy', 'wb') as f:
        np.save(f, results)
        np.save(f, lat)
        np.save(f, lon)