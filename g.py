import matplotlib.pyplot as plt
import networkx as nx
from ase.io import read
from ase.neighborlist import NeighborList, natural_cutoffs

# Load CIF file and enable periodic boundary conditions (PBC)
atoms = read("1.cif")
atoms.set_pbc(True)

cutoffs = natural_cutoffs(atoms)

nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
nl.update(atoms)

G = nx.Graph()

for i, atom in enumerate(atoms):
    G.add_node(i, element=atom.symbol, coords=atom.position)

for i in range(len(atoms)):
    indices, offsets = nl.get_neighbors(i)
    for j, offset in zip(indices, offsets):
        distance = atoms.get_distance(i, j, mic=True)  # mic=True for minimum image convention (PBC-aware)
        G.add_edge(i, j, weight=distance)

pos = {i: atoms[i].position[:2] for i in range(len(atoms))}  # Use 2D positions for visualization

plt.figure(figsize=(8, 8))

nx.draw(G, pos, with_labels=True, labels={i: G.nodes[i]['element'] for i in G.nodes()},
        node_size=500, node_color='skyblue', font_size=8, font_weight='bold')

edge_labels = nx.get_edge_attributes(G, 'weight')
nx.draw_networkx_edge_labels(G, pos, edge_labels={(i, j): f"{d:.2f} Å" for (i, j), d in edge_labels.items()},
                             font_color='red', font_size=8)

plt.title("Graph Representation of CIF File with PBC (Atoms and Bond Distances)")
plt.show()
