import ase.io
from ase.build import make_supercell
from ase.neighborlist import NeighborList, natural_cutoffs
import networkx as nx
import numpy as np

# Load the CIF file
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

P = np.eye(3) * 3
supercell = make_supercell(atoms, P)

supercell_cell = supercell.get_cell()
positions = supercell.get_scaled_positions()

central_indices = []
for i, pos in enumerate(positions):
    if all(1 / 3 <= p < 2 / 3 for p in pos):  # Central unit cell in a 3x3x3 supercell
        central_indices.append(i)

cutoffs = natural_cutoffs(supercell)
neigh_list = NeighborList(cutoffs, self_interaction=False, bothways=True)
neigh_list.update(supercell)

G = nx.Graph()
for i in range(len(supercell)):
    indices, offsets = neigh_list.get_neighbors(i)
    for idx, offset in zip(indices, offsets):
        G.add_edge(i, idx)

subgraph_list = []
visited = set()
for idx in central_indices:
    if idx not in visited:
        component = nx.node_connected_component(G, idx)
        visited.update(component)  # Mark all nodes in the component as visited
        subgraph_atoms = supercell[list(component)]
        subgraph_list.append(subgraph_atoms)
        ase.io.write(f"extracted_subgraph_{idx}.xyz", subgraph_atoms)

print(f"Extracted {len(subgraph_list)} subgraphs for central atoms.")
