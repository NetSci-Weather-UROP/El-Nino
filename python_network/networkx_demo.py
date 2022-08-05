import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network

G_krt = nx.karate_club_graph()

nt = Network('800px', '800px')
nt.from_nx(G_krt)
nt.toggle_physics(True)
nt.show_buttons(filter_=['physics'])
nt.show('nx.html')