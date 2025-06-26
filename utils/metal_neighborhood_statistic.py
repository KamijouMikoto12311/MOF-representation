import numpy as np
from sklearn.cluster import DBSCAN
from utils.load_symbols import load_symbols


def calc_metal_nonmetal_distance(supercell, metal_indices, neigh_list, clustering=True, eps=0.08, min_samples=5):
    """
    Calculate distances between metal elements and their neighboring non-metal elements,
    and compute the mean, variance, and optionally perform DBSCAN clustering on these distances.

    Parameters:
        supercell: Atoms object representing the supercell (e.g., from ASE).
        metal_indices: List of indices of metal atoms in the supercell.
        neigh_list: NeighborList object for identifying neighbors.
        clustering (bool): If True, perform DBSCAN clustering on the distances.
        eps (float): The maximum distance between two points for them to be considered as in the same neighborhood (for DBSCAN).
        min_samples (int): The number of points required to form a cluster (for DBSCAN).

    Returns:
        dict: A dictionary where keys are metal elements and values are dictionaries
              containing non-metal elements and their distances, means, variances, and optionally clusters.
              Format:
              {
                  metal_element: {
                      non_metal_element: {
                          "distances": [list of distances],
                          "mean": mean_distance,
                          "variance": variance_of_distances,
                          "clusters": [list of cluster labels] (optional, if clustering=True)
                      }
                  }
              }
    """
    element_symbols = load_symbols()
    non_metal_symbols = element_symbols["non_metal_symbols"]

    distances = {}

    # Step 1: Calculate distances between metal atoms and non-metal neighbors
    for metal_index in metal_indices:
        metal_element = supercell[metal_index].symbol
        if metal_element not in distances:
            distances[metal_element] = {}

        indices, offsets = neigh_list.get_neighbors(metal_index)

        for neighbor_index in indices:
            neighbor_index = int(neighbor_index)
            neighbor_element = supercell[neighbor_index].symbol

            if neighbor_element not in non_metal_symbols:
                continue

            if neighbor_element not in distances[metal_element]:
                distances[metal_element][neighbor_element] = {
                    "indices": [],
                    "distances": [],
                    "mean": 0,
                    "variance": 0,
                    "clusters": [],
                }

            distance = supercell.get_distance(metal_index, neighbor_index, mic=True)
            distances[metal_element][neighbor_element]["indices"].append(neighbor_index)
            distances[metal_element][neighbor_element]["distances"].append(distance)

    # Step 2: Calculate the mean and variance of distances, and optionally apply DBSCAN clustering
    for metal_element, non_metals in distances.items():
        for non_metal_element, dist_info in non_metals.items():
            dist_list = dist_info["distances"]

            mean_distance = np.mean(dist_list)
            variance_distance = np.var(dist_list)

            distances[metal_element][non_metal_element]["mean"] = mean_distance
            distances[metal_element][non_metal_element]["variance"] = variance_distance

            if clustering:
                distances_reshaped = np.array(dist_list).reshape(-1, 1)

                dbscan = DBSCAN(eps=eps, min_samples=min_samples)
                clusters = dbscan.fit_predict(distances_reshaped)

                distances[metal_element][non_metal_element]["clusters"] = clusters.tolist()

    return distances


def cluster_distances(distances):
    """
    Process clustered distances to identify the mean distances for all clusters,
    their corresponding cluster labels, and retrieve the indices belonging to each cluster.

    Parameters:
        distances: Dictionary output from the `calc_metal_nonmetal_distance` function.

    Returns:
        result: A dictionary with the same structure as distances, but including
                the mean distance, cluster label, and indices for all clusters.
                Format:
                {
                    metal_element: {
                        non_metal_element: {
                            "clusters_info": [
                                {
                                    "cluster_label": cluster_number,
                                    "mean_distance": mean_distance_of_cluster,
                                    "indices": [indices belonging to this cluster]
                                },
                                ...
                            ]
                        }
                    }
                }
    """
    result = {}

    for metal_element, non_metals in distances.items():
        result[metal_element] = {}
        for non_metal_element, dist_info in non_metals.items():
            clusters = dist_info.get("clusters", [])
            dist_list = dist_info["distances"]
            indices = dist_info["indices"]

            clusters_info = []

            for cluster_label in set(clusters):
                if cluster_label == -1:  # Ignore noise cluster
                    continue
                cluster_distances = [dist for dist, cluster in zip(dist_list, clusters) if cluster == cluster_label]
                mean_distance = np.mean(cluster_distances)

                cluster_atom_indices = [index for index, cluster in zip(indices, clusters) if cluster == cluster_label]

                clusters_info.append(
                    {"cluster_label": cluster_label, "mean_distance": mean_distance, "indices": cluster_atom_indices}
                )

            result[metal_element][non_metal_element] = {"clusters_info": clusters_info}

    return result


def get_non_max_mean_indices(result):
    """
    Get all indices of clusters that do not have the largest mean distance.

    Parameters:
        result: Output from the `process_distances` function, structured as:
                {
                    metal_element: {
                        non_metal_element: {
                            "clusters_info": [
                                {
                                    "cluster_label": cluster_number,
                                    "mean_distance": mean_distance_of_cluster,
                                    "indices": [indices belonging to this cluster]
                                },
                                ...
                            ]
                        }
                    }
                }

    Returns:
        non_max_indices: Dictionary with the same structure as `result`,
                         but containing the indices of all clusters
                         that do not have the largest mean distance.
                         Format:
                         {
                             metal_element: {
                                 non_metal_element: [indices not in the max mean cluster]
                             }
                         }
    """
    non_max_indices = {}

    for metal_element, non_metals in result.items():
        non_max_indices[metal_element] = {}
        for non_metal_element, data in non_metals.items():
            clusters_info = data.get("clusters_info", [])

            if not clusters_info:
                # No clusters, skip processing
                non_max_indices[metal_element][non_metal_element] = []
                continue

            # Find the cluster with the largest mean distance
            max_mean_cluster = max(clusters_info, key=lambda x: x["mean_distance"])
            max_mean_cluster_label = max_mean_cluster["cluster_label"]

            # Collect indices of all other clusters
            excluded_indices = []
            for cluster in clusters_info:
                if cluster["cluster_label"] != max_mean_cluster_label:
                    excluded_indices.extend(cluster["indices"])

            non_max_indices[metal_element][non_metal_element] = excluded_indices

    return non_max_indices
