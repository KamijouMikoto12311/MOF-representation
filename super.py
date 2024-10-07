import ase.io
from ase.build import make_supercell
from ase.neighborlist import NeighborList, natural_cutoffs
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import os

from utils.visualize_molecule_graph import visualize
from utils.LineGraph import create_line_graph_with_angles, visualize_line_graph

atoms = ase.io.read("str_m2_o2_o11_pcu_sym.25.cif")

metal_symbols = [
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Mn",
    "Cr",
    "V",
    "Ti",
    "Sc",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
]
non_metal_indices = [
    i for i, atom in enumerate(atoms) if atom.symbol not in metal_symbols
]
atoms = atoms[non_metal_indices]

P = np.eye(3) * 3  # Transformation matrix for 3x3x3 supercell
supercell = make_supercell(atoms, P)

positions = supercell.get_scaled_positions()

# Identify atoms that are in the central unit cell
central_indices = []
for i, pos in enumerate(positions):
    if all(1 / 3 <= p < 2 / 3 for p in pos):  # Central unit cell in a 3x3x3 supercell
        central_indices.append(i)

cutoffs = natural_cutoffs(supercell)
neigh_list = NeighborList(cutoffs, self_interaction=False, bothways=True)
neigh_list.update(supercell)

# * Create a graph(G) for supercell
G = nx.Graph()
for i in range(len(supercell)):
    G.add_node(i, element=supercell[i].symbol, position=supercell.get_positions()[i])

    indices, offsets = neigh_list.get_neighbors(i)
    for j, offset in zip(indices, offsets):
        if supercell[j].symbol not in metal_symbols:
            bond_length = supercell.get_distance(i, j, mic=True)
            G.add_edge(i, j, bond_length=bond_length)

output_dir = "ligands_xyz"
os.makedirs(output_dir, exist_ok=True)
subgraph_list = []
visited = set()
for idx in central_indices:
    if idx not in visited:
        component = nx.node_connected_component(G, idx)
        visited.update(component)  # Mark all nodes in the component as visited
        subgraph_atoms = supercell[list(component)]
        subgraph_list.append(subgraph_atoms)
        ase.io.write(os.path.join(output_dir, f"extracted_subgraph_{idx}.xyz"), subgraph_atoms)
        
        subgraph = G.subgraph(component).copy()
        L = create_line_graph_with_angles(subgraph, supercell)
        visualize(subgraph, G, idx)
        visualize_line_graph(L)
        plt.show()

print(f"Extracted {len(subgraph_list)} subgraphs for central atoms.")
