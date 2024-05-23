import numpy as np
import itertools
from Node import Node

"""Implements the UPGMA algorithm."""
class UPGMA:
    def __init__(self, dist_matrix: np.ndarray, taxa: list):
        """ Parameters:
        dist_matrix : numpy array, distance matrix of species
        taxa : list of int or str to identify taxa
        """
        self.distances = dist_matrix
        self.taxa = taxa
        self.build_tree(self.distances, self.taxa)

    def build_tree(self, dist_matrix: np.ndarray, taxa: list) -> Node:
        """Parameters:
        dist_matrix : np.ndarray of pairwise distances
        taxa : list of taxa id. Elements of lists have to be unique
        Returns the root node for constructed tree.
        """
        # Individual taxa node
        nodes = list(Node(taxon) for taxon in taxa)
        # Row/Column id to node
        rc_to_node = dict([i, j] for i, j in enumerate(nodes))
        # Taxa to row/column id
        taxa_to_rc = dict([i, j] for j, i in enumerate(taxa))
        # Copy existing distance matrix
        work_matrix = dist_matrix
        # Diaganol elements set to infinity
        work_matrix = np.array(work_matrix, dtype=float)
        np.fill_diagonal(work_matrix, np.inf)
        while len(nodes) > 1:
            # Find (row, column) of least distance
            least_id = np.unravel_index(work_matrix.argmin(), work_matrix.shape, "C")
            least_dist = work_matrix[least_id[0], least_id[1]]
            # Set node to least distance
            node1, node2 = rc_to_node[least_id[0]], rc_to_node[least_id[1]]
            # Add OTU with children node1/node2
            new_node = Node(node2, node1)
            nodes.append(new_node)
            node1.uh = least_dist / 2 - node1.dh
            node2.uh = least_dist / 2 - node2.dh
            new_node.dh = least_dist / 2
            nodes.remove(node1)
            nodes.remove(node2)
            # Update matrix
            work_matrix = self.update_distance(dist_matrix, nodes, taxa_to_rc)
            # Update row/col id to node
            rc_to_node = dict([i, j] for i, j in enumerate(nodes))
        # Set tree to root
        self.tree = nodes[0]

    def update_distance(
        self, dist_matrix: np.ndarray, nodes: list, taxa_to_rc: dict
    ) -> np.ndarray:
        """Parameters:
        dist_matrix : np.ndarray of pairwise distances for all taxa
        nodes : list of updated nodes
        taxa_to_rc : dict for taxa -> row/col id
        Returns np.ndarray of pairwise distances for updated nodes
        """
        # Node to row/col id
        node_to_rc = dict([i, j] for j, i in enumerate(nodes))
        rc = len(nodes)
        new_dist_matrix = np.zeros((rc, rc), dtype=float)
        for node1 in nodes:
            row = node_to_rc[node1]
            for node2 in nodes:
                node_pairs = list(itertools.product(node1.leaves(), node2.leaves()))
                col = node_to_rc[node2]
                new_dist_matrix[row, col] = sum(
                    dist_matrix[taxa_to_rc[i], taxa_to_rc[j]] for i, j in node_pairs
                ) / len(node_pairs)
        np.fill_diagonal(new_dist_matrix, np.inf)
        return new_dist_matrix

