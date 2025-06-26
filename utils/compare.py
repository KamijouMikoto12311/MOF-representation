from collections import Counter


def compare_graphs(graph1, graph2):
    elements_graph1 = [data["element"] for _, data in graph1.nodes(data=True)]
    elements_graph2 = [data["element"] for _, data in graph2.nodes(data=True)]

    count_graph1 = Counter(elements_graph1)
    count_graph2 = Counter(elements_graph2)

    return count_graph1 == count_graph2


def remove_duplicate(graph_list):
    """
    Remove duplicate NetworkX graphs from a list.
    Graphs are considered duplicates if they have same number of atoms for every element.
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
