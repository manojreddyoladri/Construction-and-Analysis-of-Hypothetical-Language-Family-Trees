"""Represents a node in the CMM tree"""
class Node:
    """Initializes a node"""

    def __init__(self, left=None, right=None, pass_down=[], up_height=0.0, down_height=0.0):
        """left : default = none, taxon label
        right : default = none, taxon label
        pass_down: default = Empty list, Contains the list of mutations the node passes down
        """
        self.left = left
        self.right = right
        self.pass_down = pass_down
        self.uh = up_height
        self.dh = down_height
    """Returns a list of taxa under the node"""

    def leaves(self) -> list:

        if self is None:
            return []
        if self.right is None:
            return [self.left]
        leaves = self.left.leaves() + self.right.leaves()
        return leaves

    """Returns the number of original taxa under the node"""

    def taxa_len(self) -> int:
        return sum(1 for taxa in self.leaves())

    """Returns a readable representation of the node"""

    def node_view(self) -> str:
        return "-".join(self.leaves())