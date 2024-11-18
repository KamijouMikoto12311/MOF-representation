import os
import shutil
import warnings
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from ase.io import read, write
from ase.build import make_supercell
from ase.neighborlist import NeighborList, natural_cutoffs
from concurrent.futures import ProcessPoolExecutor
from utils.visualize_molecule_graph import visualize
from utils.LineGraph import (
    create_line_graph_with_angles,
    create_line_graph_with_dihedrals,
    visualize_line_graph_with_angles,
    visualize_line_graph_with_dihedrals,
)
from utils.compare import remove_duplicate

warnings.filterwarnings("ignore")

METAL_SYMBOLS = [
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

INPUT_DIR = "cifs"
OUTPUT_DIR = "ligands_xyz"
PROCESSED_DIR = "processed_cifs"
os.makedirs(INPUT_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(PROCESSED_DIR, exist_ok=True)


def process_file(file_path):
    """Process a single file to extract molecular structures and generate visualizations."""
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    print(f"Processing >>> {file_path}")

    atoms = read(file_path)
    shutil.move(file_path, PROCESSED_DIR)

    non_metal_indices = [i for i, atom in enumerate(atoms) if atom.symbol not in METAL_SYMBOLS]
    atoms = atoms[non_metal_indices]

    P = np.eye(3) * 3
    supercell = make_supercell(atoms, P)
    positions = supercell.get_scaled_positions()

    cutoffs = natural_cutoffs(supercell)
    neigh_list = NeighborList(cutoffs, self_interaction=False, bothways=True)
    neigh_list.update(supercell)

    G = nx.Graph()
    for i in range(len(supercell)):
        G.add_node(i, element=supercell[i].symbol, position=supercell.get_positions()[i])
        indices, offsets = neigh_list.get_neighbors(i)
        for j, offset in zip(indices, offsets):
            if supercell[j].symbol not in METAL_SYMBOLS:
                bond_length = supercell.get_distance(i, j, mic=True)
                G.add_edge(i, j, bond_length=bond_length)

    central_indices = [i for i, pos in enumerate(positions) if all(1 / 3 <= p < 2 / 3 for p in pos)]

    visited = set()
    subgraph_list = []
    for idx in central_indices:
        if idx not in visited:
            component = nx.node_connected_component(G, idx)
            visited.update(component)
            subgraph = G.subgraph(component).copy()
            subgraph_list.append(subgraph)
    subgraph_list = remove_duplicate(subgraph_list)

    this_output_dir = os.path.join(OUTPUT_DIR, base_name)
    os.makedirs(this_output_dir, exist_ok=True)

    for idx, subgraph in enumerate(subgraph_list):
        this_molecule_output_dir = os.path.join(this_output_dir, f"molecule_{idx+1}")
        os.makedirs(this_molecule_output_dir, exist_ok=True)

        atom_idx = subgraph.nodes
        molecule = supercell[atom_idx]
        write(os.path.join(this_molecule_output_dir, f"molecule_{idx+1}.xyz"), molecule)

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


def run():
    file_paths = [os.path.join(INPUT_DIR, filename) for filename in os.listdir(INPUT_DIR) if os.path.isfile(os.path.join(INPUT_DIR, filename))]

    with ProcessPoolExecutor() as executor:
        executor.map(process_file, file_paths)


if __name__ == "__main__":
    run()
