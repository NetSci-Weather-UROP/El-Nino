import netCDF4 as nc
import numpy as np
import cupy as cp
# from scipy.ndimage import shift
# from cupyx.scipy.ndimage import shift
import matplotlib.pyplot as plt
from utils_cupy import *

with open('temp_data_1948_2021.npy', 'rb') as f:
    T = np.load(f)
    lat = np.load(f)
    lon = np.load(f)

T = cp.asarray(T)
lat = cp.asarray(lat)
lon = cp.asarray(lon)

def bulk_compute(years, T):
    print("Initialising dimensions and checking for errors...")
    C, T_in, T_out = year_series(T, lat, lon, 1972)  # test run

    month_days = [0, 31, 28, 31, 30, 31, 
                  30, 31, 31, 30, 31, 30]
    results = np.empty([2077, 
                        12, 2, cp.shape(T_out)[0]])

    print("The empty result array has shape: ",
           cp.shape(results))

    for shift_index in range(12):
        T = cp.roll(T, -month_days[shift_index], axis = 0)

        for year in years:
            print("Computing year:", year, "; shift:", shift_index)
            C, T_in, T_out = year_series(T, lat, lon, year)

            in_N = cp.sum((C[:,:,1] < 151), axis = 0)
            in_C = cp.sum(C[:,:,0], axis = 0)
            results[year, shift_index, 0, :] = cp.asnumpy(in_C)
            results[year, shift_index, 1, :] = cp.asnumpy(in_N)
    
    return results

years_to_compute = np.arange(1958, 1974)

results = bulk_compute(years_to_compute, T)
lat = cp.asnumpy(lat)
lon = cp.asnumpy(lon)

with open('bulk_temp.npy', 'wb') as f:
        np.save(f, results)
        np.save(f, lat)
        np.save(f, lon)