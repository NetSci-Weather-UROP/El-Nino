import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

## Taken from python implementation

#matplotlib.use('GTK3Agg')

min=-30
max=30

ds = nc.Dataset(f"generated_data.h5")
in_C = np.array(ds["in_C"])
lon = np.array(ds["lon"])
lat = np.array(ds["lat"])
lon, lat = np.meshgrid(lon, lat)

fig = plt.figure()
map = Basemap(projection='mill',lon_0=178.75)
map.pcolormesh(lon, lat, in_C, cmap="jet", latlon=True, vmin=min, vmax=max)
map.drawcoastlines()
map.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
map.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
cb = map.colorbar()
plt.show()