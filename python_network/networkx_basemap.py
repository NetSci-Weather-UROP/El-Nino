"""
Builds a NetworkX graph based on a given year's
data and attempt to implement visualisation using
Basemap.

Working directory needs to be ./el-nino/python_network

- el-nino
    - python_network
        - [networkx_cartopy.py]
    - [temp_data_1948_2021.npy]
    - [utils.py]
    - [adj.npy]*
    
*: Generated from Julia:
   https://imperiallondon-my.sharepoint.com
   /:f:/g/personal/amp320_ic_ac_uk/Erxb873d
   vJ1HlvebfsAn8T4BuEkWB6zTH2Ke8POorl5wjg?e
   =og5A6V
"""

from venv import create
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib
import sys
import inspect
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


def adjacency(C, tolerance = 1):
    """
    Outputs a relavent adjacency matrix for constructing
    graphs.

    Parameters:
    C:          3-D NumPy array containing a particular year's
                in-link information - we get it from util.py
                year_series.
    tolerance:  Determines how links are discarded according
                to the weights. See adjacency_filter from
                utils.py for more information.
    
    Returns:
    A_ij:       Weighted adjacency matrix (the el-nino basin 
                nodes first, followed by the remaining nodes).
    node_sum:   Useful.

    Note:       The construction of the adjacency matrix
                assumes there are only links between the 
                in-nodes and the out-nodes.
    """

    # get m x n matrix of weights
    C_link = adjacency_filter(C, tolerance)  # see utils
    print("Dimension of weight matrix:", np.shape(C_link))

    # get the total number of nodes (locations)
    node_sum = np.shape(C_link)[0] + np.shape(C_link)[1]
    print("node_sum:", node_sum)

    # construct A_ij
    A_ij = np.zeros([node_sum, node_sum])
    A_ij[0 : np.shape(C_link)[0],
        np.shape(C_link)[0] : node_sum] = C_link
    A_ij[np.shape(C_link)[0] : node_sum,
         0 : np.shape(C_link)[0]] = C_link.transpose()
    print("Dimension of A_ij adjacency matrix",
        np.shape(A_ij))
    
    return A_ij, node_sum


def get_pos(map, T_in, T_out, lon, lat, node_sum):
    """
    Return point_pos: dictionary of point locations to be
    plotted over a Basemap map (defined in make_graph) 
    and its unscaled original (lon + lat)
    
    Note: the positions are scaled according to the
    projection of the map.
    """
    
    lats = []
    lons = []
    for point in T_in:
        coord1, coord2 = point[-2].astype(int), point[-1].astype(int)
        lons.append(lon[coord1])
        lats.append(lat[coord2])
    for point in T_out:
        coord1, coord2 = point[-2].astype(int), point[-1].astype(int)
        lons.append(lon[coord1])
        lats.append(lat[coord2])
    
    map_x, map_y = map(lons, lats)

    point_pos = {}
    point_pos_unscaled = {}
    for i in range(node_sum):
        point_pos.update({i : (map_x[i], map_y[i])})
        point_pos_unscaled.update({i : (lons[i], lats[i])})
    
    return point_pos, point_pos_unscaled


def make_graph(year, tolerance, every_pair = False, 
               save_to_files = False):
    """
    Makes graph.

    Parameters:
    year:       The year for which the network is to be made
    tolerance:  Determines how links are discarded according
                to the weights. See adjacency_filter from
                utils.py for more information.
    every_pair: Whether we use the Python-generated data
                of only links between the in and out nodes,
                or the Julia-generated data of cross-corr.
                calculated for every pair of nodes in the
                year 1972. Enabling every_pair will lock
                the year to 1972 since it is the only
                currently available data.
    
    Returns:
    G:          The graph object itself (NetworkX).
    A_ij:       The raw weighted adjacency matrix that made
                the graph G.
    point_pos
    _unscaled:  Dictionary of point locations not scaled to 
                the Basemap.
    point_pos:  Dictionary of point locations scaled to the
                Basemap map for plotting.
    color_map:  A list to colour-code the nodes (see docstring
                above code for issues).
    """

    # set up directory to read data from parent directory
    basepath = path.dirname(__file__)
    filepath1 = path.abspath(path.join(basepath, "..",
                                       "temp_data_1948_2021.npy"))
    filepath2 = path.abspath(path.join(basepath, "..",
                                       "adj.npy"))

    # use the compiled npy data generated by nino.py
    with open(filepath1, 'rb') as f:
        T = np.load(f)
        lat = np.load(f)
        lon = np.load(f)
    
    if not every_pair:
        # get year data and calculate adjacency matrix A_ij
        print("Computing data for year:", year)
        C, T_in, T_out = year_series(T, lat, lon, year)
        print("Dimension of C matrix:", np.shape(C))
        A_ij, node_sum = adjacency(C, tolerance)
    else:
        # A_ij for every pair
        print("Computing data for year 1972",  
              "(this cannot be changed currently due",
              "to data availability (see nino.jl).")
        with open(filepath2, 'rb') as f:
            adj = np.load(f)
        print("Dimension of adj matrix:", np.shape(adj))
        node_sum = np.shape(adj)[0] * np.shape(adj)[1]
        A_ij = adj.reshape((node_sum, node_sum))
        
        # set all self-cross-corr to 0
        np.fill_diagonal(A_ij, 0)  
        print("Dimension of A_ij adjacency matrix:",
              np.shape(A_ij))

    # construct graph with weighted edges from A_ij
    G = nx.from_numpy_matrix(A_ij, parallel_edges = False)
    print("Number of edges in graph: ", G.number_of_edges())

    # initialise map
    map = Basemap(projection = 'robin',
                  lon_0 = -178.25,
                  lat_0 = 0,
                  lat_ts= 0,
                  resolution='i',  # intermediate resolution
    )

    # get node positions to plot over basemap
    point_pos, point_pos_unscaled = get_pos(map, T_in, T_out, 
                                            lon, lat, node_sum)
    
    """
    Constructs a colour map & a node size map for the in and out nodes.
    It would be a lot more useful if this were a dictionary, not a list. 
    There are ways around it though.
    """
    color_map = []
    size_map = []
    for node in G:
        if node < 57:
            color_map.append('indigo')
            size_map.append(8)
        else: 
            color_map.append('grey')
            size_map.append(0)

    # comment out when appropriate (subgraph of only 1 in-node)
    # G = G.subgraph(np.arange(56, node_sum))
    # color_map = color_map[56:node_sum]
    # size_map = size_map[56:node_sum]

    # get arrays of edges and weights separately
    edges, weights = zip(*nx.get_edge_attributes(G, "weight").items())
    edges = np.array(edges)
    weights = np.array(weights)

    # matplotlib figure
    fig, ax = plt.subplots(
        1, 1, figsize = (200, 200))

    # set up edge colouring
    cmap = plt.cm.coolwarm  # use appropriate matplotlib colour-maps
    vmin = - abs(max(weights, key = abs))
    vmax = - vmin
    print("vmin, vmax used for cmap:", vmin, vmax)
    sm = plt.cm.ScalarMappable(cmap = cmap, 
                               norm = plt.Normalize(vmin=vmin, 
                                                    vmax=vmax))
    sm.set_array([])  # set colour bar

    # draw graph over map projection
    nx.draw_networkx(G, ax = ax,
                     with_labels = False,
                     font_size = 16,
                     alpha = 0.9,
                     width = 0.5,
                     node_size = size_map,  
                     node_color = color_map,               
                     pos = point_pos,
                     edge_color = weights,
                     edge_cmap = cmap,
                     vmin = vmin,
                     vmax = vmax,
    )
    
    # misc. settings
    map.drawcoastlines()
    map.drawcountries()
    map.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
    map.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
    cbar = plt.colorbar(sm)  # draw colour bar

    # show plot and save to file
    plt.show()
    
    if save_to_files:
        print("Saving to svg...")
        fig.savefig('test.svg')
        print("Saving to png...")
        fig.savefig('test.png', dpi = 1200) 
    print("DONE!")

    return G, A_ij, point_pos, point_pos_unscaled, color_map


make_graph(1972, 2, every_pair = True, save_to_files = True)