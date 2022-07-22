import netCDF4 as nc
import numpy as np
from scipy.ndimage import shift
import matplotlib.pyplot as plt
from utils import *

with open('temp_data_1948_2021.npy', 'rb') as f:
    T = np.load(f)
    lat = np.load(f)
    lon = np.load(f)

def bulk_compute(years):
    C, T_in, T_out = year_series(T, lat, lon, 1972)  # test run

    month_days = [0, 31, 28, 31, 30, 31, 
                  30, 31, 31, 30, 31, 30]
    results = np.empty([365 * np.shape(years)[0], 
                        12, 2, np.shape(T_out)[0]])

    print(np.shape(results))

    for shift in range(12):
        T = shift(T, month_days[shift], cval = np.NaN)

        for year in years:
            print("Computing year:", year, "; shift:", shift)
            C, T_in, T_out = year_series(T, lat, lon, year)

            in_N = np.sum((C[:,:,1] < 151), axis = 0)
            in_C = np.sum(C[:,:,0], axis = 0)

            results[]

bulk_compute([1972])