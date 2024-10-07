import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from ase.neighborlist import NeighborList


def calculate_angle(vector1, vector2):
    unit_vector1 = vector1 / np.linalg.norm(vector1)
    unit_vector2 = vector2 / np.linalg.norm(vector2)
    dot_product = np.dot(unit_vector1, unit_vector2)
    angle = np.arccos(
        np.clip(dot_product, -1.0, 1.0)
    )  # Clamp to handle numerical precision issues
    return np.degrees(angle)


def create_line_graph_with_angles(G, atoms):
    L = nx.Graph()

    for edge in G.edges(data=True):
        bond = (edge[0], edge[1])
        bond = tuple(sorted(bond))  # * When adding edges, the order is sorted
        print(bond)
        L.add_node(bond, bond_length=edge[2]["weight"])

    for node in G.nodes():
        neighbors = list(G.neighbors(node))
        if len(neighbors) < 2:
            continue  # Skip if there are less than two neighbors (no angle can be formed)

        for i in range(len(neighbors)):
            for j in range(i + 1, len(neighbors)):
                bond1 = (node, neighbors[i])
                bond2 = (node, neighbors[j])

                bond1 = tuple(sorted(bond1))
                bond2 = tuple(sorted(bond2))

                pos1 = atoms[neighbors[i]].position
                pos2 = atoms[neighbors[j]].position
                pos_center = atoms[node].position

                vector1 = pos1 - pos_center  # Vector A-B
                vector2 = pos2 - pos_center  # Vector B-C

                angle = calculate_angle(vector1, vector2)

                L.add_edge(bond1, bond2, bond_angle=angle)

    # for edge in L.edges(data=True):
    #     print(
    #         f"Bonds {edge[0]} and {edge[1]} form an angle of {edge[2]['bond_angle']} degrees."
    #     )

    return L


def visualize_line_graph(L):
    pos = nx.spring_layout(L, k=1.5, iterations=50, scale=2)

    plt.figure(figsize=(12, 10))

    # Draw the nodes of the line graph (representing bonds)
    nx.draw_networkx_nodes(L, pos, node_size=500, node_color="skyblue", alpha=0.8)

    # Draw the edges (representing angles between bonds)
    nx.draw_networkx_edges(L, pos, edge_color="gray", width=2)

    # Draw the labels for the nodes (bonds between atoms)
    node_labels = {node: f"Bond {node}" for node in L.nodes()}
    nx.draw_networkx_labels(
        L, pos, labels=node_labels, font_size=10, font_color="black"
    )

    # Draw the labels for the edges (bond angles)
    edge_labels = {(u, v): f"{d['bond_angle']:.1f}°" for u, v, d in L.edges(data=True)}
    nx.draw_networkx_edge_labels(L, pos, edge_labels=edge_labels, font_size=9)

    plt.title(
        "Line Graph Representation of Molecule\n(Nodes: Bonds, Edges: Bond Angles)",
        fontsize=14,
    )
    plt.axis("off")  # Turn off axis
    plt.show()
