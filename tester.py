import matplotlib.pyplot as plt
import numpy as np
from utils import *

start_year = 1968
end_year = 1978

T, lat, lon = read_air_data(start_year)

# remove possible gap day
if not start_year % 4:
    T = T[:-1,:,:]

for year in range(start_year, end_year):

    # append yearly data
    T = np.append(T, read_air_data(year, latlon=False), axis=0)

    # remove possible gap days
    if not year % 4:
        T = T[:-1,:,:]

T0 = np.copy(T) # actual temperature measurements

T1 = avoid_seasonality(T) # seasonality removed measurements

tau_max=200
year = 1972
start_date = find_start_date(1, 7) + (year-start_year)*365
end_date = start_date + 565

T_in, T_out = separate_regions(
    T1[start_date-200:end_date,:,:], lat, lon, 190, 240, -5, 5
) # the in and out regions

n, m = np.shape(T_in)[0], np.shape(T_out)[0]

days = 365

# temporary array holding corrcoefs for all offsets
temp = np.empty([n,m,2*tau_max+1])
t_in = T_in[:,tau_max:tau_max+days]

# iterate through all offsets
for tau in range(-tau_max, tau_max+1):
    if not tau % 50:
        print("Tau:", tau)
    # set offset window for outside nodes
    t_out = T_out[:,tau_max+tau:tau_max+days+tau]

    # append corrcoefs for specific offset
    temp[:,:,tau] = pearson_coeffs(t_in, t_out)

C = np.empty([n,m,4])

# theta_i,j
C[:,:,1] = np.argmax(np.abs(temp), axis=2) - tau_max
# C(theta_i,j)
C[:,:,0] = maximod(temp, axis=2) * (C[:,:,1] < 151) * (0 <= C[:,:,1])
# mean over all theta
C[:,:,2] = np.mean(temp, axis=2)
# standard deviation over  all theta
C[:,:,3] = np.std(temp, axis=2)