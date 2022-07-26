"""
Python implementation of reproducing the paper results.
"""
import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from utils import *

def anim(T, lat, lon):
    frames = [] # for storing the generated images
    fig = plt.figure()
    lon, lat = np.meshgrid(lon, lat)
    
    map = Basemap(projection='mill',lon_0=178.75)
    
    plt.title('Temperature measurements over a year')

    for i in range(np.shape(T)[0]):
        frames.append(
            [map.pcolormesh(
                lon, lat, T[i], latlon=True, animated=True, cmap = 'seismic',
                vmin=-8, vmax=8
            )]
        )

    ani = animation.ArtistAnimation(fig, frames, interval=50,
                                    blit=True, repeat_delay=1000)

    map.drawcoastlines()
    map.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
    map.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
    cb = map.colorbar()
    cb.ax.set_ylabel('Temperature')
    plt.show()
    ani.save(f"cross_correlation_standard_temps.mkv")

    return


def plot_data(lon, lat, data, min=None, max=None, name="plot", showplots=False):

    lon, lat = np.meshgrid(lon, lat)
    fig = plt.figure()
    plt.title(f"{name}")
    map = Basemap(projection='robin',lon_0=-178.75)
    map.pcolormesh(lon, lat, data, cmap="jet", latlon=True, vmin=min, vmax=max)
    map.drawcoastlines()
    map.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
    map.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
    cb = map.colorbar()
    plt.savefig(f"./CNW-plots/{name}.png")
    
    if showplots:
        plt.show()

    return


def run(years, showplots=False):
    for year in years:
        print("Computing year", year,"...")
        C, T_in, T_out = year_series(T, lat, lon, year)

        in_degree = np.sum((0<=C[:,:,1])*(C[:,:,1] < 151), axis=0)


        # plotting
        C_plot = np.empty([71,144])
        N_plot = np.empty([71,144])
        W_plot = np.empty([71,144])
        C_plot[T_out[:,-1].astype(int), T_out[:,-2].astype(int)] = np.sum(C[:,:,0], axis=0)
        N_plot[T_out[:,-1].astype(int), T_out[:,-2].astype(int)] = in_degree
        W_plot[T_out[:,-1].astype(int), T_out[:,-2].astype(int)] = np.sum(
            (C[:,:,0]-C[:,:,2])/C[:,:,3], axis=0
        )

        plot_data(lon, lat, C_plot, min=-15, max=15, name=f"{year} in C", showplots=showplots)
        plot_data(lon, lat, N_plot, min=0, max=57, name=f"{year} in N", showplots=showplots)
        plot_data(lon, lat, W_plot, min=-100, max=100, name=f"{year} in W", showplots=showplots)
    
    return


# T, lat, lon = get_data(1948,2021)

### UNCOMMENT FOR FASTER RERUNS ON FIRST LOCAL RUN
### note that this creates a file approx 1 GB

# with open('temp_data_1948_2021.npy', 'wb') as f:
#     np.save(f, T)
#     np.save(f, lat)
#     np.save(f, lon)

with open('temp_data_1948_2021.npy', 'rb') as f:
    T = np.load(f)
    lat = np.load(f)
    lon = np.load(f)

run([1959, 1972], showplots=True)

### total runtime: ~16s 
### (~0.5s for loading and ~15.5s for running over one period window)