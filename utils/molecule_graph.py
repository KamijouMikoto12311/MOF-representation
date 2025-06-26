import networkx as nx
import matplotlib.pyplot as plt
from utils.load_symbols import load_symbols


def create_molcule_graph(supercell, neigh_list):
    """
    Creates a NetworkX graph representation of a supercell.

    Parameters:
        supercell : ASE Atoms object
            The supercell containing atom positions and elements.
        neigh_list : NeighborList object
            Object to get neighbors of atoms in the supercell.

    Returns:
        nx.Graph: A NetworkX graph where nodes represent atoms and edges represent bonds.
    """

    element_symbols = load_symbols()
    metal_symbols = set(element_symbols["metal_symbols"])
    G = nx.Graph()

    for i in range(len(supercell)):
        G.add_node(i, element=supercell[i].symbol, position=supercell.get_positions()[i])

        indices, offsets = neigh_list.get_neighbors(i)
        for j, offset in zip(indices, offsets):
            j = int(j)
            if supercell[j].symbol not in metal_symbols:
                bond_length = supercell.get_distance(i, j, mic=True)
                G.add_edge(i, j, bond_length=bond_length)

    return G


def calculate_distance_pbc(i, j, atoms):
    return atoms.get_distance(i, j, mic=True)  # mic=True ensures PBC is considered


def visualize(subgraph, G, idx):
    pos = {i: G.nodes[i]["position"][:2] for i in subgraph.nodes}

    plt.figure(1, figsize=(8, 8))
    nx.draw(
        subgraph,
        pos=pos,
        with_labels=True,
        labels={i: (subgraph.nodes[i]["element"], i) for i in subgraph.nodes()},
        node_size=500,
        node_color="lightgreen",
        font_size=8,
        font_weight="bold",
    )
    edge_labels = nx.get_edge_attributes(subgraph, "bond_length")
    try:
        nx.draw_networkx_edge_labels(
            subgraph,
            pos=pos,
            edge_labels={(i, j): f"{d:.2f} Å" for (i, j), d in edge_labels.items()},
            font_color="blue",
            font_size=8,
            connectionstyle="arc3,rad=0",
        )
    except ValueError:
        print("Too many values to unpack (expected 3)")
    plt.title(f"Subgraph {idx + 1}: No Metals")
    # plt.show()

    return
