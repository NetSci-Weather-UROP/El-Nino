"""
NetworkX -> Gephi

Working directory needs to be ./el-nino/python_network

- el-nino
    - python_network
        - [networkx_cartopy.py]
    - [temp_data_1948_2021.npy]
    - [utils.py]
    
Use plugins GeoLayout in Gephi.
"""

from venv import create
import networkx as nx
import numpy as np
import sys
import inspect
import copy
from networkx_basemap import *
from os import path  # wacky method - crutches...

"""
This is a very whacky way to import from utils.py 
in the directory above without us using a package:
https://stackoverflow.com/questions/714063/importing-modules-
from-parent-folder
"""
currentdir = path.dirname(path.abspath(
                          inspect.getfile(
                          inspect.currentframe())))
parentdir = path.dirname(currentdir)
sys.path.insert(0, parentdir)
from utils import *  # import utils

def split_graph(G):
    '''
    Split the original weighted network into two
    according to the edge weights (positive/negative).
    The resulting graphs will both have only positive
    weights.
    '''
    
    # get edge weights (+/-)
    edge_weights = nx.get_edge_attributes(G,'weight')
    
    # create copies (doesn't take much time)
    G_pos_weights = copy.deepcopy(G)
    G_neg_weights = copy.deepcopy(G)
    
    # filter edges
    G_pos_weights.remove_edges_from((e for e, w in 
                                     edge_weights.items() if w < 0))
    G_neg_weights.remove_edges_from((e for e, w in 
                                     edge_weights.items() if w > 0))
    
    # set negative weights to be positive
    for u, v, d in G_neg_weights.edges(data = True):
        d['weight'] = - d['weight']
    
    return G_pos_weights, G_neg_weights

# make relavent graph (see networkx_basemap.py)
year = 1972
tolerance = 2
G, A_ij, point_pos, point_pos_unscaled, color_map = make_graph(year, 
                                                               tolerance)
# DEFUNCT
# gephi_loc = np.empty([np.shape(A_ij)[0], 3])
# for i in range(np.shape(A_ij)[0]):
#    gephi_loc[i] = [i, point_pos_unscaled[i][0],
#                    point_pos_unscaled[i][1]]
# np.savetxt("gephi_loc.csv", gephi_loc, delimiter = ",", fmt = '%f', 
# header = 'Id, longitude, latitude', comments='')

# store attributes to nodes
for i in range(np.shape(A_ij)[0]):
    G.nodes[i]['longitude'] = point_pos_unscaled[i][0]
    G.nodes[i]['latitude'] = point_pos_unscaled[i][1]
    
G_pos_weights, G_neg_weights = split_graph(G)

# save to GEXF
nx.write_gexf(G, path = f"{year}_tol{tolerance}.gexf")
nx.write_gexf(G_pos_weights, path 
              = f"{year}_tol{tolerance}_pos.gexf")
nx.write_gexf(G_neg_weights, path 
              = f"{year}_tol{tolerance}_neg.gexf")