import numpy as np
import networkx as nx
from math import cos, sin, isclose

with open("adj", "rb") as f:
    adj = np.load(f).transpose()

with open("lat", "rb") as f:
    lat = np.load(f)

with open("lon", "rb") as f:
    lon = np.load(f)

l_lat = len(lat)
l_lon = len(lon)

G = nx.DiGraph()

r = 15
r_lat = range(r, l_lat - r, r+1)
r_lon = range(r, l_lon - r, r+1)

print(f"{l_lat} {l_lon}")
for i in r_lat:
    for j in r_lon:
        p_lat = np.mean(lat[(i-r):(i+r+1)])
        p_lon = np.mean(lon[(j-r):(j+r+1)])
        G.add_node(f"{p_lat}:{p_lon}",
            lat = p_lat,
            lon = p_lon
        )

for i in r_lat:
    print(f"lat: {i}")
    for j in r_lon:
        for i_alt in r_lat:
            for j_alt in r_lon:
                t = np.mean(adj[
                    (i-r):(i+r+1),
                    (j-r):(j+r+1),
                    (i_alt-r):(i_alt+r+1),
                    (j_alt-r):(j_alt+r+1)
                ])
                if not(isclose(t, 0.)):
                    p_lat = np.mean(lat[(i-r):(i+r+1)])
                    p_lon = np.mean(lon[(j-r):(j+r+1)])
                    
                    p_lat_alt = np.mean(lat[(i_alt-r):(i_alt+r+1)])
                    p_lon_alt = np.mean(lon[(j_alt-r):(j_alt+r+1)])

                    G.add_edge(
                        f"{p_lat}:{p_lon}",
                        f"{p_lat_alt}:{p_lon_alt}",
                        weight = t
                    )

nx.write_gexf(G, "weather.gexf.gz")
