"""
Builds a NetworkX graph based on a given year's
data and attempt to implement visualisation using
Cartopy.

Working directory needs to be ./el-nino/python_network

- el-nino
    - python_network
        - [networkx_cartopy.py]
    - [temp_data_1948_2021.npy]
    - [utils.py]

Install Cartopy w. Conda:
conda install -c conda-forge cartopy
"""

import networkx as nx
import numpy as np
import cartopy.crs as ccrs
from cartopy.io import shapereader as shpreader
from pyvis.network import Network
import matplotlib.pyplot as plt
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


def make_graph(year = 1972, tolerance = 1):
    """
    Makes graph.

    Parameters:
    year:       The year for which the network is to be made
    tolerance:  Determines how links are discarded according
                to the weights. See adjacency_filter from
                utils.py for more information.
    
    Returns:
    G:          The graph object itself (NetworkX).
    A_ij:       The raw weighted adjacency matrix that made
                the graph G.
    point_pos:  Dictionary of point locations to be
                plotted over a cartopy map later (see docstring
                above code for issues).
    color_map:  A list to colour-code the nodes (see docstring
                above code for issues).
    """

    # set up directory to read data from parent directory
    basepath = path.dirname(__file__)
    filepath = path.abspath(path.join(basepath, "..",
                                      "temp_data_1948_2021.npy"))

    # use the compiled npy data generated by nino.py
    with open(filepath, 'rb') as f:
        T = np.load(f)
        lat = np.load(f)
        lon = np.load(f)
    
    # get year data and calculate adjacency matrix A_ij
    year = 1972
    print("Computing data for year:", year)
    C, T_in, T_out = year_series(T, lat, lon, year)
    print("Dimension of C matrix:", np.shape(C))
    A_ij, node_sum = adjacency(C, tolerance)

    # construct graph with weighted edges from A_ij
    G = nx.from_numpy_matrix(A_ij, parallel_edges = False)
    print("Number of edges in graph: ", G.number_of_edges())

    """
    Construct point_pos: dictionary of point locations to be
    plotted over a cartopy map later.

    Note: the order of lon and lat needs to be swapped here.
    The nodes are identified with a number, e.g.
    53 : [lon, lat] in point_pos refers to the node number 53
    (a node in the el-nino basin as it is within the first
    57 nodes) and its location.

    Issues: it just does not really work yet. The points would
    be placed correctly but the networkx plotting of edges
    won't follow with the new node positions...
    """
    point_pos = {}
    i = 0
    for point in T_in:
        coord1, coord2 = point[-2].astype(int), point[-1].astype(int)
        # I am attempting to use the following method to line up the
        # points on the map (see docstring for issue)
        point_pos.update({i:[(lon[coord1] + 178.25) % 360, lat[coord2]]})
        i += 1
    for point in T_out:
        coord1, coord2 = point[-2].astype(int), point[-1].astype(int)
        point_pos.update({i:[(lon[coord1] + 178.25) % 360, lat[coord2]]})
        i += 1

    """
    Constructs a colour map to colour-code the in and out nodes.
    Currently not used in this script - it would be a lot more useful
    if this were a dictionary, not a list. There are ways
    around it though.
    """
    color_map = []
    for node in G:
        if node < 57:
            color_map.append('red')
        else: 
            color_map.append('blue')

    # comment out when appropriate (subgraph of only 1 in-node)
    G = G.subgraph(np.arange(56, node_sum))

    # draw background map with cartopy
    crs = ccrs.PlateCarree(central_longitude = 178.25)
    fig, ax = plt.subplots(
        1, 1, figsize=(200, 200),
        subplot_kw = dict(projection = crs))
    ax.coastlines()
    ax.set_global()

    # draw graph with networkx over cartopy
    nx.draw_networkx(G, ax = ax,
                     font_size = 16,
                     alpha = 0.8,
                     width = 0.5,
                     node_size = 8,
                     with_labels = False,
                     pos = point_pos
    )

    plt.show()
    fig.savefig('test.png')

    return G, A_ij, point_pos, color_map


make_graph(1972, 2)