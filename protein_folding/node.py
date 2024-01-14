import protein_folding.definitions as definitions
from .vector import Vec3D

from typing import Optional


_bond_values = {
    frozenset({'H', 'H'}): -1,
    frozenset({'C', 'C'}): -5,
    frozenset({'C', 'H'}): -1,
}

_delta_pos_from_direction = {
    definitions.UP: Vec3D(0, 1, 0),
    definitions.DOWN: Vec3D(0, -1, 0),
    definitions.LEFT: Vec3D(-1, 0, 0),
    definitions.RIGHT: Vec3D(1, 0, 0),
    definitions.FORWARD: Vec3D(0, 0, 1),
    definitions.BACKWARD: Vec3D(0, 0, -1),
}


class NotNeighbourError(Exception):
    pass


class Node:
    counter = 0

    def __init__(self, letter, x, y, z, direction: int | None, prev_node: Optional['Node'] = None):
        self.id = Node.counter
        Node.counter += 1

        self.letter = letter
        self.pos = Vec3D(x, y, z)

        self.next: Node | None = None
        self.prev: Node | None = prev_node

        self.direction_from_previous = direction

        self.ghost: bool = False

    def __eq__(self, other):
        return self.id == other.id

    def __str__(self):
        return f"id: {self.id}"

    def __repr__(self):
        return (f'Node(id={self.id}, letter="{self.letter}", (x,y,z)=({self.x}, {self.y}, {self.z}), '
                f'direction={self.direction_from_previous})')

    @classmethod
    def from_previous(cls, c, direction: int, prev: 'Node'):
        x, y, z = calc_position_from_direction(direction, prev)

        return cls(c, x, y, z, direction, prev_node=prev)

    @property
    def x(self):
        return self.pos.x

    @property
    def y(self):
        return self.pos.y

    @property
    def z(self):
        return self.pos.z

    def is_neightbour(self, other: 'Node') -> bool:
        if abs(self.id - other.id) == 1:
            return False

        delta = Vec3D.abs_diff(self.pos, other.pos)

        return delta.sum_components() == 1

    def touch_direction(self, other: 'Node'):
        if not self.is_neightbour(other):
            raise NotNeighbourError(f"Nodes {str(self)} and {str(other)} are not neighbours")

        delta = Vec3D.abs_diff(self.pos, other.pos)

        return delta.eq_components(1)

    def bond_value(self, other: 'Node'):
        return _bond_values.get(frozenset({self.letter, other.letter}), 0)

    def change_direction(self, pos_set: set[Vec3D], direction: int):
        if not self.direction_from_previous:
            raise Exception("First node in chain!")

        new_pos = calc_position_from_direction(direction, self.prev)
        delta_pos = new_pos - self.pos

        # Remove old position from the positions set
        pos_set.remove(self.pos)

        self.pos = new_pos

        self.direction_from_previous = direction

        if self.next:
            self.next.cascade_position(pos_set, delta_pos)

        # Add new position to the positions set
        pos_set.add(self.pos)

        # Set node to not-ghost
        self.ghost = False

    def cascade_position(self, pos_set: set[Vec3D], delta_pos: Vec3D):
        """
        Updates the position of current node because the previous node's
            position was updated. Calls the next node and repeates the
            process. Also removes itself from the protein's positions set
            if necessary to prevent position conflicts

        pre:
            - pos_set is a set of 3D vectors representing all not-ghosted
                nodes' positions
            - delta_pos is a 3D vector representing the change in position
                compared to the node's current position

        post:
            - the node is ghosted if it wasn't ghosted yet
            - the node's position is updated
        """
        # If node wasn't ghosted, delete its position from the protein class'
        # positions set to prevent overwriting values
        if not self.ghost:
            pos_set.remove(self.pos)

        self.pos += delta_pos

        if self.next is not None:
            self.next.cascade_position(pos_set, delta_pos)


def calc_position_from_direction(direction: int, prev: Node):
    if direction not in _delta_pos_from_direction:
        raise Exception(f"Invalid direction of {direction}")

    return prev.pos + _delta_pos_from_direction[direction]
