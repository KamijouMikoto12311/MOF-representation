import networkx as nx
import matplotlib.pyplot as plt


def calculate_distance_pbc(i, j, atoms):
    return atoms.get_distance(i, j, mic=True)  # mic=True ensures PBC is considered


def visualize(subgraph, G, idx, fc=False):
    pos = {i: G.nodes[i]["position"][:2] for i in subgraph.nodes}

    plt.figure(1, figsize=(8, 8))
    nx.draw(
        subgraph,
        pos,
        with_labels=True,
        labels={i: (subgraph.nodes[i]["element"], i) for i in subgraph.nodes()},
        node_size=500,
        node_color="lightgreen",
        font_size=8,
        font_weight="bold",
    )
    edge_labels = nx.get_edge_attributes(subgraph, "bond_length")
    nx.draw_networkx_edge_labels(
        subgraph,
        pos,
        edge_labels={(i, j): f"{d:.2f} Å" for (i, j), d in edge_labels.items()},
        font_color="blue",
        font_size=8,
    )
    plt.title(f"Subgraph {idx + 1} Before Completion: No Metals")
    # plt.show()

    return
