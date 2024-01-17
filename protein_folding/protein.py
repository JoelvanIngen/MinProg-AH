import matplotlib.pyplot as plt
from .definitions import *
from .node import Node
from .vector import *

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

        # Ensure valid-ness of sequence
        _validate_protein_letters(sequence)

        self.sequence = sequence

        # Create list of all the nodes
        self.nodes = [Node(self, 0, self.sequence[0], 0, 0, 0, direction=None)]
        for _id, c in enumerate(self.sequence[1:], start=1):
            # Initialise new node in a straight line
            self.nodes.append(Node.from_previous(self, _id, c, RIGHT, self.nodes[-1]))

        # Dict to keep track of Node positions
        self.pos_to_node: dict[Vec3D, Node] = {}
        self.collect_node_positions()

    def __len__(self) -> int:
        return len(self.nodes)

    def get_order(self):
        order = [node.direction_from_previous for node in self.nodes[1:]]
        return order

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
            node.change_direction(direction, ignore_pos_set=True)

        self.collect_node_positions()

    def calc_size_score(self):
        """
        Computes a score based on the physical dimensions of the protein. This
            value can be used as a penalty for a specific configuration.
        """
        dim = get_min_max([node.pos for node in self.nodes])
        box = dim[1] - dim[0]

        return box.len_sq()

    def calc_area(self):
        """
        Computes the area of the grid, only taking into account x- and
            y-coordinates.

        post: an integer is returned representing the area occupied by the protein
        """
        dim = get_min_max([node.pos for node in self.nodes])
        box = dim[1] - dim[0]

        return box.area()

    def calc_volume(self):
        """
        Computes the physical volume of the protein, taking into account
            three dimensions.

        post: an integer is returned representing the volume occupied by the protein
        """
        dim = get_min_max([node.pos for node in self.nodes])
        box = dim[1] - dim[0]

        return box.volume()

    def collect_node_positions(self):
        self.pos_to_node = {}
        for node in self.nodes:
            self.pos_to_node[node.pos] = node

    def get_bond_score(self) -> float:
        """
        Computes the quality of the current fold purely based on the amount
            of desirable bonds and their values

        pre:
        post:
            - a float is returned representing the bond score
        """

        quality = 0

        neighbours = self.get_all_neighbours()
        for pairing in neighbours:
            node1, node2 = pairing
            quality += node1.bond_value(node2)

        return quality

    def get_all_neighbours(self) -> list[tuple[Node, Node]]:
        neighbours: list[tuple[Node, Node]] = []
        for i, node1 in enumerate(self.nodes[:-1]):
            for node2 in self.nodes[i + 1:]:
                if node1.is_neighbour(node2):
                    neighbours.append((node1, node2))

        return neighbours

    def get_node_id(self, node: Node):
        for i in range(len(self.nodes)):
            if node == self.nodes[i]:
                return i

        raise KeyError

    def has_valid_order(self) -> bool:
        """
        Determines if every node has a unique position.

        post:
            - Returns True if every node has a unique position
            - Returns False if two nodes overlap
        """

        taken_positions = set()

        for node in self.nodes:
            if node.pos in taken_positions:
                return False

            taken_positions.add(node.pos)

        return True

    def filter_neighbours_by_nonzero_score(self, neighbours: list[tuple[Node, Node]]):
        """
        Finds and removes neighbours whose bond score is zero from the
            list of neighbours.
        """

        filtered_neighbours: list[tuple[Node, Node]] = []
        for (node1, node2) in neighbours:
            if node1.bond_value(node2):
                filtered_neighbours.append((node1, node2))

        return filtered_neighbours

    def straighten(self):
        """
        Sets protein order to be a straight line from left to right
        """
        self.set_order([RIGHT] * (len(self) - 1))

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

        prev = Vec3D(0, 0, 0)

        for n in self.nodes:
            # Put node letter at node's position
            plt.text(n.x, n.y, n.letter, size='10')

            # Draw protein line segment from previous node
            plt.plot([prev.x, n.x], [prev.y, n.y], '-', color='black', linewidth=2)

            # Save position as previous
            prev = n.pos

        # Add lines between all pairings
        neighbours = self.get_all_neighbours()
        neighbours_filtered = self.filter_neighbours_by_nonzero_score(neighbours)
        for pairing in neighbours_filtered:
            node1, node2 = pairing

            if node1.letter == 'H' and node2.letter == 'H':
                line_colour = 'green'
            elif (node1.letter == 'H' and node2.letter == 'C'
                  or node1.letter == 'C' and node2.letter == 'H'):
                line_colour = 'orange'
            elif node1.letter == 'C' and node2.letter == 'C':
                line_colour = 'red'
            else:
                # Should never trigger, for now this mostly shows nothing gets through the if-statements
                raise Exception(f"Letters {node1.letter} and {node2.letter} should not be neighbours")

            plt.plot([node1.x, node2.x], [node1.y, node2.y], '--', color=line_colour, linewidth=1)

        # Get full dimensions of protein
        dim = get_min_max(self.nodes)

        plt.xlim(dim[0].x - 1, dim[1].x + 1)
        plt.ylim(dim[0].y - 1, dim[1].y + 1)
        plt.axis('off')

        # Save image
        plt.savefig(filename)


def _validate_protein_letters(seq: str) -> None:
    """
    Checks whether a sequence is valid

    pre:
        - seq is a str of protein letters

    post:
        - InvalidSequenceError is raised if any letter in the sequence is
        invalid
    """
    for c in seq:
        if c not in _valid_protein_letters:
            raise InvalidSequenceError(f"Letter {c} in sequence {seq} is not valid")
