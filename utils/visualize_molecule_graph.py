import networkx as nx
import matplotlib.pyplot as plt


def calculate_distance_pbc(i, j, atoms):
    return atoms.get_distance(i, j, mic=True)  # mic=True ensures PBC is considered


def visualize(subgraph, G, idx, fc=False):
    """
    Visualize a molecular subgraph and optionally complete it to a fully connected graph.

    This function visualizes a molecular subgraph, represented as a NetworkX graph object,
    by displaying nodes and edges. The visualization can occur in two stages:
    1. Visualization of the original subgraph, highlighting existing connections.
    2. Optionally completing the subgraph into a fully connected graph and visualizing it.

    Parameters:
    -----------
    subgraph : networkx.Graph
        The subgraph to visualize. Nodes represent atoms, and edges represent bonds or interactions.

    G : networkx.Graph
        The parent graph from which the subgraph was extracted. Used to obtain node positions.

    idx : int
        The index of the subgraph, used for titling the plots.

    fc : bool, optional, default=False
        If True, the subgraph is completed to a fully connected graph, and the completed version
        is visualized in addition to the original subgraph.

    Visualization Details:
    ----------------------
    - Nodes are colored 'lightgreen' for the original subgraph and 'skyblue' for the fully connected version.
    - Nodes are labeled with their element type and index in the original subgraph.
    - Edges are labeled with bond distances (in Å) for both visualizations.
    - Metals are excluded from both visualizations.

    Notes:
    ------
    - The function uses a distance calculation function (`calculate_distance_pbc`) to account for
      Periodic Boundary Conditions (PBC) when completing the subgraph.
    - The visualization highlights the organic atoms and their connections in the molecular graph.

    Returns:
    --------
    None

    Example Usage:
    --------------
    >>> visualize(subgraph, G, idx=0, fc=True)
    This will visualize the given subgraph, followed by the completed fully connected version.
    """
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

    if fc:
        # * Complete the subgraph to a fully connected graph
        nodes = list(subgraph.nodes)
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                if not subgraph.has_edge(nodes[i], nodes[j]):
                    distance = calculate_distance_pbc(nodes[i], nodes[j], atoms)
                    subgraph.add_edge(nodes[i], nodes[j], weight=distance)

        plt.figure(3, figsize=(8, 8))
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
        # plt.show()

    return
