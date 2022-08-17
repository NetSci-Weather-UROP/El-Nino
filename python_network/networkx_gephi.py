"""
NetworkX -> Gephi

Working directory needs to be ./el-nino/python_network

- el-nino
    - python_network
        - [networkx_cartopy.py]
    - [temp_data_1948_2021.npy]
    - [utils.py]
"""

from venv import create
import networkx as nx
import numpy as np
import sys
import inspect
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

make_graph(1972, 2)