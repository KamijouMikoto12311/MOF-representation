import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from ase.neighborlist import NeighborList


def calculate_angle(vector1, vector2):
    unit_vector1 = vector1 / np.linalg.norm(vector1)
    unit_vector2 = vector2 / np.linalg.norm(vector2)
    dot_product = np.dot(unit_vector1, unit_vector2)
    angle = np.arccos(np.clip(dot_product, -1.0, 1.0))  # Clamp to handle numerical precision issues
    return np.degrees(angle)


def calculate_dihedral(vector1, vector2, vector3, cif_name, processing_idx):
    normal1 = np.cross(vector1, vector2)
    normal2 = np.cross(vector2, vector3)

    if np.linalg.norm(normal1) == 0 or np.linalg.norm(normal2) == 0:
        if bool(processing_idx):
            print(
                f"{cif_name} >>> molecule_{processing_idx}: Zero vector encountered in dihedral calculation. Probably because of 180° bond angle. Corresponding dihedrals shown as nan."
            )
        else:
            print(f"{cif_name} >>> Zero vector encountered in dihedral calculation. Probably because of 180° bond angle. Corresponding dihedrals shown as nan.")
        return np.nan

    normal1 /= np.linalg.norm(normal1)
    normal2 /= np.linalg.norm(normal2)

    dot_product = np.dot(normal1, normal2)
    dot_product = np.clip(dot_product, -1.0, 1.0)

    angle = np.arccos(dot_product)
    return 180 - np.degrees(angle)


def create_line_graph_with_angles(G, atoms):
    L = nx.Graph()

    for edge in G.edges(data=True):
        bond = (edge[0], edge[1])
        bond = tuple(sorted(bond))  # * When adding edges, the order is sorted
        L.add_node(bond, bond_length=edge[2]["bond_length"])

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
                vector2 = pos2 - pos_center  # Vector C-B

                angle = calculate_angle(vector1, vector2)

                L.add_edge(bond1, bond2, bond_angle=angle)

    return L


def create_line_graph_with_dihedrals(L, atoms, all=True, cif_name=None, processing_idx=0):
    """
    Creates a second-level line graph where:
        - Nodes represent angles (edges in L).
        - Edges represent dihedral angles between those angles.
    Parameters:
        L (networkx.Graph): Line graph with bond angles as attributes.
        atoms (ASE Atoms object): Atomic structure with positions.
    Returns:
        LL (networkx.Graph): Line graph of dihedrals.
    """
    LL = nx.Graph()

    for edge in L.edges(data=True):  # * Each edge of L represents an angle
        angle = (edge[0], edge[1])  # Bonds forming the angle
        angle = tuple(sorted(angle))  # Consistent ordering
        LL.add_node(angle, bond_angle=edge[2]["bond_angle"])

    for bond in L.nodes():  # * Each node of L represents a bond
        neighbors = list(L.neighbors(bond))
        if len(neighbors) < 2:
            continue

        for i in range(len(neighbors)):
            for j in range(i + 1, len(neighbors)):
                angle1 = (bond, neighbors[i])
                angle2 = (bond, neighbors[j])

                if bool(set(neighbors[i]) & set(neighbors[j])) and not all:
                    continue

                else:
                    a, b = [(a, b) for a in range(len(bond)) for b in range(len(neighbors[i])) if bond[a] == neighbors[i][b]][0]
                    _, d = [(c, d) for c in range(len(bond)) for d in range(len(neighbors[j])) if bond[c] == neighbors[j][d]][0]
                    vector_1 = atoms[neighbors[i][0 if b == 1 else 1]].position - atoms[neighbors[i][b]].position
                    vector_share = atoms[bond[0 if a == 1 else 1]].position - atoms[bond[a]].position
                    vector_2 = atoms[neighbors[j][0 if d == 1 else 1]].position - atoms[neighbors[j][d]].position

                    dihedral = calculate_dihedral(vector_1, vector_share, vector_2, cif_name, processing_idx)

                    angle1 = tuple(sorted(angle1))
                    angle2 = tuple(sorted(angle2))

                    LL.add_edge(angle1, angle2, dihedral_angle=dihedral)
    return LL


def visualize_line_graph_with_angles(L):
    pos = nx.spring_layout(L, k=1.5, iterations=50, scale=2)

    plt.figure(2, figsize=(12, 10))

    nx.draw_networkx_nodes(L, pos, node_size=500, node_color="skyblue", alpha=0.8)
    nx.draw_networkx_edges(L, pos, edge_color="gray", width=2)
    node_labels = {node: f"Bond {node}" for node in L.nodes()}
    nx.draw_networkx_labels(L, pos, labels=node_labels, font_size=10, font_color="black")
    edge_labels = {(u, v): f"{d['bond_angle']:.1f}°" for u, v, d in L.edges(data=True)}
    nx.draw_networkx_edge_labels(L, pos, edge_labels=edge_labels, font_size=9)
    plt.title(
        "Line Graph Representation of Molecule\n(Nodes: Bonds, Edges: Bond Angles)",
        fontsize=14,
    )
    plt.axis("off")


def visualize_line_graph_with_dihedrals(LL):
    pos = nx.spring_layout(LL, k=2.0, iterations=50, scale=2)

    plt.figure(3, figsize=(12, 10))

    nx.draw_networkx_nodes(LL, pos, node_size=500, node_color="lightgreen", alpha=0.8)
    nx.draw_networkx_edges(LL, pos, edge_color="purple", width=2)
    node_labels = {node: f"Angle {node}" for node in LL.nodes()}
    nx.draw_networkx_labels(LL, pos, labels=node_labels, font_size=10, font_color="black")
    edge_labels = {(u, v): f"{d['dihedral_angle']:.1f}°" for u, v, d in LL.edges(data=True)}
    nx.draw_networkx_edge_labels(LL, pos, edge_labels=edge_labels, font_size=9)

    plt.title(
        "Second-Level Line Graph (Dihedral Representation)\n(Nodes: Angles, Edges: Dihedrals)",
        fontsize=14,
    )
    plt.axis("off")
