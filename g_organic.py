import matplotlib.pyplot as plt
import networkx as nx
from ase.io import read, write
from ase.neighborlist import NeighborList, natural_cutoffs
from sklearn.manifold import MDS
import numpy as np


# Function to calculate distance considering PBC
def calculate_distance_pbc(i, j, atoms):
    return atoms.get_distance(i, j, mic=True)  # mic=True ensures PBC is considered


def save_subgraph_to_xyz(subgraph, idx):
    # Extract distance matrix based on edge weights
    nodes = list(subgraph.nodes)
    num_nodes = len(nodes)
    dist_matrix = np.zeros((num_nodes, num_nodes))

    # Fill in the distance matrix
    for i, node_i in enumerate(nodes):
        for j, node_j in enumerate(nodes):
            if i == j:
                dist_matrix[i, j] = 0.0  # Diagonal is zero (same atom)
            elif subgraph.has_edge(node_i, node_j):
                dist_matrix[i, j] = subgraph.edges[node_i, node_j]["weight"]
            else:
                raise Exception

    # Use MDS to reconstruct the 3D coordinates
    mds = MDS(n_components=3, dissimilarity="precomputed", random_state=42)
    coords = mds.fit_transform(dist_matrix)

    # Output to XYZ format
    xyz_filename = f"subgraph_{idx + 1}_mds_reconstructed.xyz"
    with open(xyz_filename, "w") as f:
        f.write(f"{num_nodes}\n")
        f.write(f"Subgraph {idx + 1} reconstructed coordinates using MDS\n")
        for i, node in enumerate(nodes):
            element = subgraph.nodes[node]["element"]
            x, y, z = coords[i]
            f.write(f"{element} {x:.6f} {y:.6f} {z:.6f}\n")

    print(f"Subgraph {idx + 1} saved to {xyz_filename}")


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

atoms = read("str_m2_o2_o11_pcu_sym.25.cif")
atoms.set_pbc(True)

# Generate natural cutoffs for neighbor detection
cutoffs = natural_cutoffs(atoms)

# Create NeighborList, with PBC enabled
nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
nl.update(atoms)

# Create an empty graph
G = nx.Graph()

for i, atom in enumerate(atoms):
    if atom.symbol not in metal_elements:
        G.add_node(i, element=atom.symbol, coords=atom.position)

# Add edges (bonds) based on neighbor list, considering PBC, but only between non-metal atoms
for i in range(len(atoms)):
    if atoms[i].symbol in metal_elements:
        continue

    indices, offsets = nl.get_neighbors(i)
    for j, offset in zip(indices, offsets):
        if atoms[j].symbol not in metal_elements:
            distance = atoms.get_distance(i, j, mic=True)
            G.add_edge(i, j, weight=distance)

# Extract subgraphs
subgraphs = [G.subgraph(c).copy() for c in nx.connected_components(G)]


# Visualize and process each subgraph
for idx, subgraph in enumerate(subgraphs):
    pos = {
        i: G.nodes[i]["coords"][:2] for i in subgraph.nodes
    }  # Use 2D positions for visualization

    # Plot the initial (incomplete) subgraph
    plt.figure(figsize=(8, 8))
    nx.draw(
        subgraph,
        pos,
        with_labels=True,
        labels={i: subgraph.nodes[i]["element"] for i in subgraph.nodes()},
        node_size=500,
        node_color="lightgreen",
        font_size=8,
        font_weight="bold",
    )
    edge_labels = nx.get_edge_attributes(subgraph, "weight")
    nx.draw_networkx_edge_labels(
        subgraph,
        pos,
        edge_labels={(i, j): f"{d:.2f} Å" for (i, j), d in edge_labels.items()},
        font_color="blue",
        font_size=8,
    )
    plt.title(f"Subgraph {idx + 1} Before Completion: No Metals")
    plt.show()

    # Complete the subgraph to a fully connected graph
    nodes = list(subgraph.nodes)
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            if not subgraph.has_edge(nodes[i], nodes[j]):
                distance = calculate_distance_pbc(nodes[i], nodes[j], atoms)
                subgraph.add_edge(nodes[i], nodes[j], weight=distance)

    save_subgraph_to_xyz(subgraph, idx)

    # Plot the completed (fully connected) subgraph
    plt.figure(figsize=(8, 8))
    nx.draw(
        subgraph,
        pos,
        with_labels=True,
        labels={i: subgraph.nodes[i]["element"] for i in subgraph.nodes()},
        node_size=500,
        node_color="skyblue",
        font_size=8,
        font_weight="bold",
    )
    edge_labels = nx.get_edge_attributes(subgraph, "weight")
    nx.draw_networkx_edge_labels(
        subgraph,
        pos,
        edge_labels={(i, j): f"{d:.2f} Å" for (i, j), d in edge_labels.items()},
        font_color="red",
        font_size=8,
    )
    plt.title(
        f"Fully Connected Subgraph {idx + 1}: Organic Atoms and Bond Distances (No Metals)"
    )
    plt.show()
