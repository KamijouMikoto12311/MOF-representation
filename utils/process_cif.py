import os
import json
import time
import shutil
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from ase.io import read, write
from ase.build import make_supercell
from ase.neighborlist import NeighborList, natural_cutoffs

from utils.molecule_graph import visualize, create_molcule_graph
from utils.line_graph import (
    create_line_graph_with_angles,
    create_line_graph_with_dihedrals,
    visualize_line_graph_with_angles,
    visualize_line_graph_with_dihedrals,
)
from utils.compare import remove_duplicate
from utils.load_symbols import load_symbols


def get_non_metals_with_only_metal_neighbors(supercell, neigh_list):
    element_symbols = load_symbols()
    metal_symbols = element_symbols["metal_symbols"]
    non_metal_symbols = element_symbols["non_metal_symbols"]

    non_metals_with_metal_neighbors = []

    for i, atom in enumerate(supercell):
        if atom.symbol in non_metal_symbols:
            neighbor_indices, _ = neigh_list.get_neighbors(i)

            if all(supercell[neighbor_idx].symbol in metal_symbols for neighbor_idx in neighbor_indices):
                non_metals_with_metal_neighbors.append(i)

    return non_metals_with_metal_neighbors


def process_cif(file_path, linegraph=False):
    """Process a single file to extract molecular structures and generate visualizations."""
    output_dir = "./data/ligands"
    processed_dir = "./data/processed_cifs"
    neglected_dir = "./data/neglected_cifs"
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(processed_dir, exist_ok=True)

    element_symbols = load_symbols()
    metal_symbols = element_symbols["metal_symbols"]
    non_metal_symbols = element_symbols["non_metal_symbols"]

    cif_name = os.path.splitext(os.path.basename(file_path))[0]
    atoms = read(file_path)
    atoms_num = len(atoms)

    if atoms_num > 1000:
        print(f"Neglecting >>> {cif_name}, more than 1000 atoms in cif.")
        os.makedirs(neglected_dir, exist_ok=True)
        shutil.move(file_path, neglected_dir)
        return

    else:
        print(f"Processing >>> {cif_name}.")

        supercell = make_supercell(atoms, np.eye(3) * 3)
        metal_indices = [i for i, atom in enumerate(supercell) if atom.symbol in metal_symbols]
        cutoffs = natural_cutoffs(supercell)
        neigh_list = NeighborList(cutoffs, skin=0.1, self_interaction=False, bothways=True)
        neigh_list.update(supercell)

        non_metals_with_only_metal_neighbors = get_non_metals_with_only_metal_neighbors(supercell, neigh_list)
        remove_indices = non_metals_with_only_metal_neighbors + metal_indices

        supercell = supercell[[i for i in range(len(supercell)) if i not in remove_indices]]
        cutoffs = natural_cutoffs(supercell)
        neigh_list = NeighborList(cutoffs, skin=0.1, self_interaction=False, bothways=True)
        neigh_list.update(supercell)

        G = create_molcule_graph(supercell, neigh_list)

        central_indices = [
            int(i) for i, pos in enumerate(supercell.get_scaled_positions()) if all(1 / 3 <= p < 2 / 3 for p in pos)
        ]

        visited = set()
        subgraph_list = []
        for idx in central_indices:
            if idx not in visited:
                component = nx.node_connected_component(G, idx)
                visited.update(component)
                subgraph = G.subgraph(component).copy()
                subgraph_list.append(subgraph)
        subgraph_list = remove_duplicate(subgraph_list)

        this_output_dir = os.path.join(output_dir, cif_name)
        os.makedirs(this_output_dir, exist_ok=True)

        for idx, subgraph in enumerate(subgraph_list):
            this_molecule_output_dir = os.path.join(this_output_dir, f"molecule_{idx+1}")
            os.makedirs(this_molecule_output_dir, exist_ok=True)

            atom_idx = subgraph.nodes
            molecule = supercell[atom_idx]
            write(os.path.join(this_molecule_output_dir, f"molecule_{idx+1}.xyz"), molecule)  # * Write the xyz file
            visualize(subgraph, G, idx)
            plt.savefig(
                os.path.join(this_molecule_output_dir, "mol_graph.png"), format="png"
            )  # * Save the molecule figure
            plt.close()

            if linegraph:
                L = create_line_graph_with_angles(subgraph, supercell)
                LL = create_line_graph_with_dihedrals(L, supercell, all=True, cif_name=cif_name, processing_idx=idx)

                visualize_line_graph_with_angles(L)
                plt.savefig(os.path.join(this_molecule_output_dir, "line_graph.png"), format="png")
                plt.close()
                visualize_line_graph_with_dihedrals(LL)
                plt.savefig(os.path.join(this_molecule_output_dir, "line_line_graph.png"), format="png")
                plt.close()

        shutil.move(file_path, processed_dir)
