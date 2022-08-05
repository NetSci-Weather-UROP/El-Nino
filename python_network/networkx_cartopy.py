import networkx as nx
import numpy as np
import cartopy.crs as ccrs
from cartopy.io import shapereader as shpreader
import matplotlib.pyplot as plt

# testing
edges = np.array([['A', 'B'],
                 ['B', 'C'],
                 ['C', 'D'],
                 ['C', 'A']])
pos = {
       'A': [150, 5], 
       'B': [110, 15], 
       'C': [140, -25], 
       'D': [160, -30]
}

G = nx.from_edgelist(edges)
sg = next(G.subgraph(c) for c in nx.connected_components(G))

# draw background map
crs = ccrs.PlateCarree(central_longitude = 178.75)
fig, ax = plt.subplots(
    1, 1, figsize=(10.24, 7.68),
    subplot_kw=dict(projection=crs))
ax.coastlines()
ax.set_global()

nx.draw_networkx(sg, ax=ax,
                 font_size=16,
                 alpha=1,
                 width=2,
                 pos=pos,
                 cmap=plt.cm.autumn)
                 
plt.show()