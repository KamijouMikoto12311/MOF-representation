import matplotlib.pyplot as plt
import networkx as nx
from ase.io import read, write
from ase.neighborlist import NeighborList, natural_cutoffs
from sklearn.manifold import MDS
import numpy as np

from utils.LineGraph import create_line_graph_with_angles, visualize_line_graph
from utils.visualize_molecule_graph import visualize


# Define a list of metal elements (this can be extended as needed)
metal_elements = [
    "Cu",
    "Zn",
    "Fe",
    "Co",
    "Ni",
    "Mn",
    "Ag",
    "Al",
    "Mg",
    "Ti",
    "Zr",
]

atoms = read("str_m7_o4_o24_bcu_sym.62.cif")
atoms.set_pbc(True)

cutoffs = natural_cutoffs(atoms)

nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
nl.update(atoms)

G = nx.Graph()

for i, atom in enumerate(atoms):
    if atom.symbol not in metal_elements:
        G.add_node(i, element=atom.symbol, position=atom.position)

# Add edges (bonds) based on neighbor list, considering PBC, but only between non-metal atoms
for i in range(len(atoms)):
    if atoms[i].symbol in metal_elements:
        continue

    indices, offsets = nl.get_neighbors(i)
    for j, offset in zip(indices, offsets):
        if atoms[j].symbol not in metal_elements:
            distance = atoms.get_distance(i, j, mic=True)
            G.add_edge(i, j, bond_length=distance)

subgraphs = [G.subgraph(c).copy() for c in nx.connected_components(G)]

# Visualize and process each subgraph
for idx, subgraph in enumerate(subgraphs):
    L = create_line_graph_with_angles(subgraph, atoms)
    visualize(subgraph, G, idx=idx, fc=True)
    visualize_line_graph(L)
