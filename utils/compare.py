import networkx as nx
from networkx.algorithms import isomorphism


def compare_graphs(graph1, graph2):
    """
    Compare two NetworkX graphs to check if they are isomorphic.
    The node attribute 'element' is compared as well.
    """
    node_match = isomorphism.categorical_node_match("element", None)

    return nx.is_isomorphic(graph1, graph2, node_match=node_match)


def remove_duplicate(graph_list):
    """
    Remove duplicate NetworkX graphs from a list.
    Graphs are considered duplicates if they are isomorphic (same structure and node 'element' attributes).
    Returns a new list of unique graphs.
    """
    unique_graphs = []

    for graph in graph_list:
        is_duplicate = False
        for unique_graph in unique_graphs:
            if compare_graphs(graph, unique_graph):
                is_duplicate = True
                break

        if not is_duplicate:
            unique_graphs.append(graph)

    return unique_graphs
