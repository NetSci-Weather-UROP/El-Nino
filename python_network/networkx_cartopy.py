import networkx as nx
import numpy as np
import cartopy.crs as ccrs
from cartopy.io import shapereader as shpreader
import matplotlib.pyplot as plt

# draw background map
crs = ccrs.PlateCarree(central_longitude = 178.75)
fig, ax = plt.subplots(
    1, 1, figsize=(10.24, 7.68),
    subplot_kw=dict(projection=crs))
ax.coastlines()
ax.set_global()

plt.show()