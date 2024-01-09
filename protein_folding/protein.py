import numpy as np
import matplotlib.pyplot as plt
from protein_folding.definitions import *
from protein_folding.node import Node

_valid_protein_letters = {'H', 'P', 'C'}


class InvalidSequenceError(Exception):
    pass


class Protein:

    def __init__(self, sequence: str) -> None:
        """
        representation of a simple protein consisting of polar or hydrophobic
        amino acids in a chain.

        pre:
            - sequence should be string containing only letter H and P
        post:
            - protein is initialised in straight line ordering
        """
        # Ensure sequence only consists of valid letters
        if not all(c in _valid_protein_letters for c in sequence):
            raise InvalidSequenceError
        self.sequence = sequence

        # Create list of all the nodes
        self.nodes = [Node(self.sequence[0], 0, 0, 0, direction=None)]

        for c in self.sequence[1:]:
            # Initalise new node in a straight line
            self.nodes.append(Node.from_previous(c, RIGHT, self.nodes[-1]))

    def set_order(self, order: list) -> None:
        """
        sets order (shape) of protein to provided order, and recalculates
        coordinates of all amino acids under new order. Orders are formatted
        as a list of 1, -1, 2, -2 in some order, corresponding to left, right,
        up and down respectively, ie: [1,1,-2,-1,2] (len = len(seq) - 1)

        pre:
            - order is a list formatted as described
        post:
            - attributes for acid coordinates and order are updated
        """

        assert len(self.nodes) == len(order) + 1, \
            f"Wrong order size, got {len(order)} but expected {(len(self.nodes) - 1)}"

        for node, direction in zip(self.nodes[1:], order):
            node.change_direction(direction)

    def calc_volume(self, area_only=False):
        coords = [(node.x, node.y, node.z) for node in self.nodes]

        x = [coord[0] for coord in coords]
        y = [coord[1] for coord in coords]
        z = [coord[2] for coord in coords]

        dx = max(x) - min(x)
        dy = max(y) - min(y)
        dz = max(z) - min(z)

        if area_only:
            return dx * dy

        return dx * dy * dz

    def get_order_quality(self) -> float:
        """
        quantises the quality of the current fold, or order. Depends on bonds
        and compactness (?)

        pre:
        post:
            - a float is returned representing the fold quality
        """

        quality = 0
        for i, node1 in enumerate(self.nodes[:-1]):
            for node2 in self.nodes[i + 1:]:
                if node1.is_neightbour(node2):
                    quality += node1.bond_value(node2)

        return quality

    def plot(self, filename="./unnamed_protein.png") -> None:
        """
        save current protein configuration as image.
        pre:
            - filename is a string referring to a valid storage location and
            ending in a valid image format
        post:
            - order is saved as image under filename
        """
        fig = plt.figure()
        x, y = 0, 0
        for i, acid in enumerate(self.sequence):
            x_new, y_new = self.acid_coords[i]
            # place character
            plt.text(x_new, y_new, acid, size='10')
            # place line
            plt.plot([x, x_new], [y, y_new], '--', color='black', linewidth=1)
            x, y = x_new, y_new

        # plot settings
        plt.xlim(-len(self.sequence), len(self.sequence))
        plt.ylim(-len(self.sequence), len(self.sequence))
        plt.axis('off')
        # save image
        plt.savefig(filename)


if __name__ == "__main__":
    seq = "PHPPHPHH"
    order = [2, -1, 2, 1, 1, -2, -2]
    protein = Protein(seq)
    protein.set_order(order)
    print(f"Order quality: {protein.get_order_quality()}")  # should be -2.0
    protein.plot()
