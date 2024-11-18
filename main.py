import ase.io
from ase.build import make_supercell
from ase.neighborlist import NeighborList, natural_cutoffs
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
import shutil

from utils.visualize_molecule_graph import visualize
from utils.line_graph import (
    create_line_graph_with_angles,
    create_line_graph_with_dihedrals,
    visualize_line_graph_with_angles,
    visualize_line_graph_with_dihedrals,
)
from utils.compare import remove_duplicate

warnings.filterwarnings("ignore")

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

input_dir = "cifs"
output_dir = "ligands"
processed_dir = "processed_cifs"
os.makedirs(input_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)
os.makedirs(processed_dir, exist_ok=True)

for filename in os.listdir(input_dir):
    file_path = os.path.join(input_dir, filename)
    base_name = os.path.splitext(filename)[0]
    print(f"Processing >>> {file_path}")

    atoms = ase.io.read(file_path)
    shutil.move(file_path, processed_dir)

    non_metal_indices = [i for i, atom in enumerate(atoms) if atom.symbol not in metal_symbols]
    atoms = atoms[non_metal_indices]

    P = np.eye(3) * 3  # Transformation matrix for 3x3x3 supercell
    supercell = make_supercell(atoms, P)
    positions = supercell.get_scaled_positions()

    cutoffs = natural_cutoffs(supercell)
    neigh_list = NeighborList(cutoffs, self_interaction=False, bothways=True)
    neigh_list.update(supercell)

    # * Create a graph(G) for supercell. Note that the index of node is the same as the atom index in supercell!!
    G = nx.Graph()
    for i in range(len(supercell)):
        G.add_node(i, element=supercell[i].symbol, position=supercell.get_positions()[i])

        indices, offsets = neigh_list.get_neighbors(i)
        for j, offset in zip(indices, offsets):
            if supercell[j].symbol not in metal_symbols:
                bond_length = supercell.get_distance(i, j, mic=True)
                G.add_edge(i, j, bond_length=bond_length)

    molecule_list = []
    subgraph_list = []
    visited = set()
    # * Identify atoms that are in the central unit cell
    central_indices = []
    for i, pos in enumerate(positions):
        if all(1 / 3 <= p < 2 / 3 for p in pos):
            central_indices.append(i)

    for idx in central_indices:
        if idx not in visited:
            component = nx.node_connected_component(G, idx)
            visited.update(component)  # Mark all nodes in the component as visited

            subgraph = G.subgraph(component).copy()
            subgraph_list.append(subgraph)

    subgraph_list = remove_duplicate(subgraph_list)

    this_output_dir = os.path.join(output_dir, base_name)
    os.makedirs(this_output_dir, exist_ok=True)

    for idx, subgraph in enumerate(subgraph_list):
        this_molecule_output_dir = os.path.join(this_output_dir, f"molecule_{idx+1}")
        os.makedirs(this_molecule_output_dir, exist_ok=True)

        atom_idx = subgraph.nodes
        molecule = supercell[atom_idx]
        ase.io.write(os.path.join(this_molecule_output_dir, f"molecule_{idx+1}.xyz"), molecule)

        L = create_line_graph_with_angles(subgraph, supercell)
        LL = create_line_graph_with_dihedrals(L, supercell, all=True, cif_name=base_name, processing_idx=idx)

        visualize(subgraph, G, idx)
        plt.savefig(os.path.join(this_molecule_output_dir, "mol_graph.png"), format="png")
        plt.close()
        visualize_line_graph_with_angles(L)
        plt.savefig(os.path.join(this_molecule_output_dir, "line_graph.png"), format="png")
        plt.close()
        visualize_line_graph_with_dihedrals(LL)
        plt.savefig(os.path.join(this_molecule_output_dir, "line_line_graph.png"), format="png")
        plt.close()
