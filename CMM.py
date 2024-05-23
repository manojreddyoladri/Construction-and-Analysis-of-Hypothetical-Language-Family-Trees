import itertools
from Node import Node


"""Implements the CMM algorithm."""
class CMM:
    def __init__(self, cmm: list, taxa: list):
        """ Parameters:
        cmm : cmm matrix
        taxa : list of str to identify taxa
        """
        self.cmm = cmm
        self.taxa = taxa
        self.build_tree(self.cmm, self.taxa)

    def build_tree(self, cmm: list, taxa: list) -> Node:
        """Parameters:
        cmm: cmm matrix
        taxa : list of taxa id. Elements of lists have to be unique
        Returns the root node for constructed tree.
        """
        # Individual taxa node
        nodes = list(Node(taxon) for taxon in taxa)
        # Row/Column id to node
        rc_to_node = dict([i, j] for i, j in enumerate(nodes))
        # Taxa to row/column id
        taxa_to_rc = dict([i, j] for j, i in enumerate(taxa))
        # Copy existing cmm matrix
        work_matrix = cmm
        count = 0
        while(True):
            max_row_index = -1
            max_col_index = -1
            max_list_length = 0
            for i, row in enumerate(work_matrix):
                for j, lst in enumerate(row):
                    if i != j and lst is not None:  # Exclude diagonal elements
                        current_list_length = len(lst)
                        if current_list_length > max_list_length:
                            max_list_length = current_list_length
                            max_row_index = i
                            max_col_index = j
            if(max_row_index == -1 or len(work_matrix[max_row_index][max_col_index])==0):
                break
            count+=1
            node1, node2 = rc_to_node[max_row_index], rc_to_node[max_col_index]
            # Add OTU with children node1/node2
            new_node = Node(node2, node1)
            new_node.pass_down = work_matrix[max_row_index][max_col_index]
            node1.pass_down = [elem for elem in (node1.pass_down or []) if elem not in new_node.pass_down]
            node2.pass_down = [elem for elem in (node2.pass_down or []) if elem not in new_node.pass_down]
            nodes.append(new_node)
            nodes.remove(node1)
            nodes.remove(node2)
            # Update matrix
            work_matrix = self.update_cmm(cmm, nodes, taxa_to_rc)
            # Update row/col id to node
            rc_to_node = dict([i, j] for i, j in enumerate(nodes))
        # Set tree to root
        self.tree = nodes

    def update_cmm(self, cmm: list, nodes: list, taxa_to_rc: dict
    ) -> list:
        # Node to row/col id
        node_to_rc = dict([i, j] for j, i in enumerate(nodes))
        rc = len(nodes)
        new_cmm = [[None for i in range(rc)] for j in range(rc)]
        for node1 in nodes:
            row = node_to_rc[node1]
            for node2 in nodes:
                node_pairs = list(itertools.product(node1.leaves(), node2.leaves()))
                col = node_to_rc[node2]
                i,j = node_pairs[0]
                new_cmm[row][col] = cmm[taxa_to_rc[i]][taxa_to_rc[j]]
                for idx in range(1,len(node_pairs)):
                    i,j = node_pairs[idx][0], node_pairs[idx][1]
                    # Check if new_cmm[row][col] is not None before attempting to iterate
                    if new_cmm[row][col] is not None and cmm[taxa_to_rc[i]][taxa_to_rc[j]] is not None:
                        new_cmm[row][col] = list(set(new_cmm[row][col]) & set(cmm[taxa_to_rc[i]][taxa_to_rc[j]]))
                    else:
                        new_cmm[row][col] = cmm[taxa_to_rc[i]][taxa_to_rc[j]] if cmm[taxa_to_rc[i]][taxa_to_rc[j]] is not None else None
        return new_cmm
