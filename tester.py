import matplotlib.pyplot as plt
import numpy as np
from utils import *

start_year = 1948
end_year = 2021

# T, lat, lon = read_air_data(start_year)

# # remove possible gap day
# if not start_year % 4:
#     T = T[:-1,:,:]

# for year in range(start_year, end_year):

#     # append yearly data
#     T = np.append(T, read_air_data(year, latlon=False), axis=0)

#     # remove possible gap days
#     if not year % 4:
#         T = T[:-1,:,:]

# T0 = np.copy(T) # actual temperature measurements

# T1 = avoid_seasonality(T) # seasonality removed measurements

with open('temp_data_1948_2021.npy', 'rb') as f:
    T1 = np.load(f)
    lat = np.load(f)
    lon = np.load(f)

tau_max=200

# year = 2009
# loc = "Mozambique"
# outx, outy = 13, 44

# year = 2002
# loc = "India"
# outx, outy = 30, 24

year = 1972
loc = "Australia"
outx, outy = 57, 45

# year = 1963
# loc = "Sea of Okhotsk"
# outx, outy = 59, 14


print(loc, lon[outx], lat[outy])

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
    temp[:,:,tau_max+tau] = pearson_coeffs(t_in, t_out)

C = np.empty([n,m,4])

# theta_i,j
C[:,:,1] = np.argmax(np.abs(temp), axis=2) - tau_max
# C(theta_i,j)
C[:,:,0] = maximod(temp, axis=2) * (C[:,:,1] < 151) * (0 <= C[:,:,1])
# mean over all theta
C[:,:,2] = np.mean(temp, axis=2)
# standard deviation over  all theta
C[:,:,3] = np.std(temp, axis=2)

M = reformat_c(C, T_out)

inx, iny = 92, 36

print("in-location", lon[inx], lat[iny])

fig = plt.figure()
plt.title(f"{year} anomalies at {loc}, {outx*2.5}° lon {87.5-outy*2.5}° lat")
plt.plot(T1[start_date:start_date+365,iny,inx])
plt.plot(T1[start_date:start_date+365,outy,outx])
plt.savefig(f"./CNW-plots/{year}-{loc}-anomalies.png")
plt.show()

outind = list(
    set(
        np.where(T_out[:,-2]==outx)[0]
    ).intersection(set(np.where(T_out[:,-1]==outy)[0]))
)[0]

inind = list(
    set(
        np.where(T_in[:,-2]==inx)[0]
    ).intersection(set(np.where(T_in[:,-1]==iny)[0]))
)[0]

print(
    "out location", T_out[outind, -2]*2.5, 87.5-T_out[outind, -1]*2.5, "\n in location", T_in[inind, -2]*2.5, 87.5-T_in[inind, -1]*2.5
)

# fig = plt.figure()
# plt.title(f"{year} crosscoef progression")
# plt.ylim([-0.5,0.6])
# for i in range(57):
#     plt.plot(np.arange(-200,201),temp[i,outind,::-1],linewidth=0.4)
# # plt.savefig(f"./CNW-plots/{year}-crosscoef-progression.png")
# plt.show()

fig = plt.figure()
plt.title(f"{year} crosscoef in {loc}")
plt.ylim([-0.4,0.4])
plt.plot(np.arange(-200,201),temp[inind,outind,:])
plt.savefig(f"./CNW-plots/{year}-crosscoef-{loc}.png")
plt.show()