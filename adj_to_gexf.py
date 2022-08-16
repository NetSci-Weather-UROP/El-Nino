import numpy as np
import networkx as nx
from math import cos, sin, isclose
from matplotlib import pyplot as plt

with open("adj", "rb") as f:
    adj = np.load(f).transpose()

with open("lat", "rb") as f:
    lat = np.load(f)

with open("lon", "rb") as f:
    lon = np.load(f)

l_lat = len(lat)
l_lon = len(lon)

G = nx.Graph()

for i in range(l_lat):
    for j in range(l_lon):
        G.add_node(f"{lat[i]}:{lon[j]}", pos={
            "x" : cos(lat[i]) * cos(lon[j]),
            "y" : cos(lat[i]) * sin(lon[j]),
            "z" : sin(lat[i])
        })

for i in range(l_lat):
    print(f"lat: {i}")
    for j in range(l_lon):
        for i_alt in range(l_lat):
            for j_alt in range(l_lon):
                t = adj[i,j,i_alt,j_alt]
                if not(isclose(t, 0.)):
                    G.add_edge(
                        f"{lat[i]}:{lon[j]}",
                        f"{lat[i_alt]}:{lon[j_alt]}",
                        weight = t
                    )

nx.write_gexf(G, "weather.gexf.gz")
