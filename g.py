import matplotlib.pyplot as plt
import networkx as nx
from ase.io import read
from ase.neighborlist import NeighborList, natural_cutoffs

# Load CIF file and create graph
atoms = read("str_m2_o2_o11_pcu_sym.25.cif")
cutoffs = natural_cutoffs(atoms)
nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
nl.update(atoms)

# Create an empty graph
G = nx.Graph()

# Add atoms as nodes, labeled by their element symbols
for i, atom in enumerate(atoms):
    G.add_node(i, element=atom.symbol, coords=atom.position)

# Add edges based on neighbor list (bond detection), labeled by distance (bond length)
for i in range(len(atoms)):
    indices, offsets = nl.get_neighbors(i)
    for j in indices:
        distance = atoms.get_distance(i, j)
        G.add_edge(i, j, weight=distance)

# Prepare for visualization: Node positions (2D projection for visualization)
pos = {i: atoms[i].position[:2] for i in range(len(atoms))}  # Use 2D positions for visualization

# Draw the graph
plt.figure(figsize=(8, 8))

# Draw the nodes with labels (element symbols)
nx.draw(G, pos, with_labels=True, labels={i: G.nodes[i]['element'] for i in G.nodes()},
        node_size=500, node_color='skyblue', font_size=8, font_weight='bold')

# Draw edges (bonds) with labels (distances)
edge_labels = nx.get_edge_attributes(G, 'weight')
nx.draw_networkx_edge_labels(G, pos, edge_labels={(i, j): f"{d:.2f} Å" for (i, j), d in edge_labels.items()},
                             font_color='red', font_size=8)

plt.title("Graph Representation of CIF File (Atoms and Bond Distances)")
plt.show()
