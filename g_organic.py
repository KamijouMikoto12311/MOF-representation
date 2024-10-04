import matplotlib.pyplot as plt
import networkx as nx
from ase.io import read
from ase.neighborlist import NeighborList, natural_cutoffs
import numpy as np

# Define a list of metal elements (this can be extended as needed)
metal_elements = ["Cu", "Zn", "Fe", "Co", "Ni", "Mn", "Ag", "Al", "Mg", "Ti", "Zr"]  # Add any other metals in your structure

# Load CIF file and enable periodic boundary conditions (PBC)
atoms = read("str_m2_o2_o11_pcu_sym.25.cif")
atoms.set_pbc(True)  # Enable periodic boundary conditions (PBC)

# Generate natural cutoffs for neighbor detection
cutoffs = natural_cutoffs(atoms)

# Create NeighborList, with PBC enabled
nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
nl.update(atoms)

# Create an empty graph
G = nx.Graph()

# Add only non-metal atoms as nodes, labeled by their element symbols
for i, atom in enumerate(atoms):
    if atom.symbol not in metal_elements:  # Exclude metals
        G.add_node(i, element=atom.symbol, coords=atom.position)

# Add edges (bonds) based on neighbor list, considering PBC, but only between non-metal atoms
for i in range(len(atoms)):
    if atoms[i].symbol in metal_elements:
        continue  # Skip metal atoms

    indices, offsets = nl.get_neighbors(i)
    for j, offset in zip(indices, offsets):
        if atoms[j].symbol not in metal_elements:  # Only connect non-metal atoms
            distance = atoms.get_distance(i, j, mic=True)  # mic=True for minimum image convention (PBC-aware)
            G.add_edge(i, j, weight=distance)

# Find all connected components (subgraphs)
subgraphs = [G.subgraph(c).copy() for c in nx.connected_components(G)]

# Function to calculate Euclidean distance between two 3D points
def calculate_distance(coords1, coords2):
    return np.linalg.norm(np.array(coords1) - np.array(coords2))

# Visualize each fully connected subgraph
for idx, subgraph in enumerate(subgraphs):
    pos = {i: G.nodes[i]['coords'][:2] for i in subgraph.nodes}  # Use 2D positions for visualization
    
    # Complete the subgraph to a fully connected graph
    nodes = list(subgraph.nodes)
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            if not subgraph.has_edge(nodes[i], nodes[j]):  # If edge doesn't already exist
                # Calculate the distance between the two nodes (atoms)
                coords_i = subgraph.nodes[nodes[i]]['coords']
                coords_j = subgraph.nodes[nodes[j]]['coords']
                distance = calculate_distance(coords_i, coords_j)
                
                # Add the edge with the calculated distance as the weight
                subgraph.add_edge(nodes[i], nodes[j], weight=distance)
    
    # Plot the fully connected subgraph
    plt.figure(figsize=(8, 8))

    # Draw the nodes with labels (element symbols)
    nx.draw(subgraph, pos, with_labels=True, labels={i: subgraph.nodes[i]['element'] for i in subgraph.nodes()},
            node_size=500, node_color='skyblue', font_size=8, font_weight='bold')

    # Draw edges (bonds) with labels (distances)
    edge_labels = nx.get_edge_attributes(subgraph, 'weight')
    nx.draw_networkx_edge_labels(subgraph, pos, edge_labels={(i, j): f"{d:.2f} Å" for (i, j), d in edge_labels.items()},
                                 font_color='red', font_size=8)

    plt.title(f"Fully Connected Subgraph {idx + 1}: Organic Atoms and Bond Distances (No Metals)")
    plt.show()
