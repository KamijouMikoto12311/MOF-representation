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


def calculate_dihedral(vector1, vector2, vector3):
    """
    Calculate the dihedral angle between three vectors.
    """
    # Normal vectors of the planes
    normal1 = np.cross(vector1, vector2)
    normal2 = np.cross(vector2, vector3)

    # Check for zero vectors, which would lead to NaN
    if np.linalg.norm(normal1) == 0 or np.linalg.norm(normal2) == 0:
        print("Zero vector encountered in dihedral calculation.")
        return np.nan

    # Normalize the normals
    normal1 /= np.linalg.norm(normal1)
    normal2 /= np.linalg.norm(normal2)

    # Calculate the angle between the two normal vectors
    dot_product = np.dot(normal1, normal2)
    # Clamp the dot product to avoid numerical issues with arccos
    dot_product = np.clip(dot_product, -1.0, 1.0)

    angle = np.arccos(dot_product)

    # Determine the sign using the direction of vector2
    sign = np.dot(normal1, vector3)
    if sign < 0:
        angle = -angle

    print(f"Dihedral angle calculated: {np.degrees(angle):.2f}°")
    return np.degrees(angle)


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
                vector2 = pos2 - pos_center  # Vector B-C

                angle = calculate_angle(vector1, vector2)

                L.add_edge(bond1, bond2, bond_angle=angle)

    return L


# def create_line_graph_with_dihedrals(L, atoms):
#     LL = nx.Graph()

#     for edge in L.edges(data=True):  # * edge = angle
#         angle = (edge[0], edge[1])  # * angle = Tuple(bond1, bond2)
#         angle = tuple(sorted(angle))  # * When adding edges, the order is sorted
#         LL.add_node(angle, bond_angle=edge[2]["bond_angle"])

#     for node in L.nodes():  # * node = bond
#         neighbors = list(L.neighbors(node))
#         if len(neighbors) < 2:
#             continue  # Skip if there are less than two neighbors (no angle can be formed)

#         for i in range(len(neighbors)):
#             for j in range(i + 1, len(neighbors)):
#                 angle1 = (node, neighbors[i])
#                 angle2 = (node, neighbors[j])

#                 angle1 = tuple(sorted(angle1))
#                 angle2 = tuple(sorted(angle2))

#                 # TODO: stopped at this point

#                 pos1 = atoms[neighbors[i]].position
#                 pos2 = atoms[neighbors[j]].position
#                 pos_center = atoms[node].position

#                 vector1 = pos1 - pos_center  # Vector A-B
#                 vector2 = pos2 - pos_center  # Vector B-C

#                 angle = calculate_angle(vector1, vector2)

#                 LL.add_edge(angle1, angle2, bond_angle=angle)

#     return LL


def create_line_graph_with_dihedrals(L, atoms):
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

    # Add nodes to LL, representing angles
    for edge in L.edges(data=True):  # Each edge represents an angle
        angle = (edge[0], edge[1])  # Bonds forming the angle
        angle = tuple(sorted(angle))  # Consistent ordering
        LL.add_node(angle, bond_angle=edge[2]["bond_angle"])

    # Add edges to LL, representing dihedrals between angles
    for node in L.nodes():  # Each node represents a bond
        neighbors = list(L.neighbors(node))

        for i in range(len(neighbors)):
            for j in range(i + 1, len(neighbors)):
                # Create angle1 and angle2
                angle1 = (node, neighbors[i])
                angle2 = (node, neighbors[j])

                angle1 = tuple(sorted(angle1))
                angle2 = tuple(sorted(angle2))

                # Gather all atoms from angle1 and angle2
                all_atoms = [
                    angle1[0][0],
                    angle1[0][1],
                    angle1[1][0],
                    angle2[0][0],
                    angle2[0][1],
                    angle2[1][0],
                ]

                # Remove duplicates by converting to a set, then back to a list
                unique_atoms = list(dict.fromkeys(all_atoms))

                # Check if we have exactly four unique atoms
                if len(unique_atoms) != 4:
                    continue  # Skip if we don't have four distinct atoms

                # Assign the four unique atoms to A, B, C, D
                pos_a = atoms[unique_atoms[0]].position
                pos_b = atoms[unique_atoms[1]].position
                pos_c = atoms[unique_atoms[2]].position
                pos_d = atoms[unique_atoms[3]].position

                # Calculate vectors
                vector_ab = pos_b - pos_a  # A -> B
                vector_bc = pos_c - pos_b  # B -> C
                vector_cd = pos_d - pos_c  # C -> D

                # Check for zero vectors before calculating dihedral
                if (
                    np.linalg.norm(vector_ab) == 0
                    or np.linalg.norm(vector_bc) == 0
                    or np.linalg.norm(vector_cd) == 0
                ):
                    print("Zero vector encountered in dihedral calculation. Skipping.")
                    continue

                # Calculate dihedral angle
                dihedral_angle = calculate_dihedral(vector_ab, vector_bc, vector_cd)

                # Add edge to the second-level line graph LL
                LL.add_edge(angle1, angle2, dihedral_angle=dihedral_angle)

    return LL


def visualize_line_graph_with_angles(L):
    pos = nx.spring_layout(L, k=1.5, iterations=50, scale=2)

    plt.figure(2, figsize=(12, 10))

    nx.draw_networkx_nodes(L, pos, node_size=500, node_color="skyblue", alpha=0.8)
    nx.draw_networkx_edges(L, pos, edge_color="gray", width=2)
    node_labels = {node: f"Bond {node}" for node in L.nodes()}
    nx.draw_networkx_labels(
        L, pos, labels=node_labels, font_size=10, font_color="black"
    )
    edge_labels = {(u, v): f"{d['bond_angle']:.1f}°" for u, v, d in L.edges(data=True)}
    nx.draw_networkx_edge_labels(L, pos, edge_labels=edge_labels, font_size=9)
    plt.title(
        "Line Graph Representation of Molecule\n(Nodes: Bonds, Edges: Bond Angles)",
        fontsize=14,
    )
    plt.axis("off")


def visualize_line_graph_with_dihedrals(LL):
    pos = nx.spring_layout(LL, k=1.5, iterations=50, scale=2)

    plt.figure(3, figsize=(12, 10))

    nx.draw_networkx_nodes(LL, pos, node_size=500, node_color="lightgreen", alpha=0.8)
    nx.draw_networkx_edges(LL, pos, edge_color="purple", width=2)
    node_labels = {node: f"Angle {node}" for node in LL.nodes()}
    nx.draw_networkx_labels(
        LL, pos, labels=node_labels, font_size=10, font_color="black"
    )
    edge_labels = {
        (u, v): f"{d['dihedral_angle']:.1f}°" for u, v, d in LL.edges(data=True)
    }
    nx.draw_networkx_edge_labels(LL, pos, edge_labels=edge_labels, font_size=9)

    plt.title(
        "Second-Level Line Graph (Dihedral Representation)\n(Nodes: Angles, Edges: Dihedrals)",
        fontsize=14,
    )
    plt.axis("off")
