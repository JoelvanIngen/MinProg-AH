import definitions

from typing import Optional


_bond_values = {
    frozenset({'H', 'H'}): 1,
    frozenset({'C', 'C'}): 5,
    frozenset({'C', 'H'}): 1,
}


class NotNeighbourError(Exception):
    pass


class Node:
    def __init__(self, letter, x, y, z, direction: int | None, prev_node: Optional['Node'] = None):
        self.letter = letter
        self.x: int = x
        self.y: int = y
        self.z: int = z

        self.next: Node | None = None
        self.prev: Node | None = prev_node

        self.direction_from_previous = direction

        self._neighbours: list[Node]

    @classmethod
    def from_previous(cls, c, direction: int, prev: 'Node'):
        x = prev.x
        y = prev.y
        z = prev.z

        match direction:
            case definitions.UP:
                y += 1
            case definitions.DOWN:
                y -= 1
            case definitions.LEFT:
                x -= 1
            case definitions.RIGHT:
                x += 1

        return cls(c, x, y, z, direction, prev_node=prev)

    def is_neightbour(self, other: 'Node') -> bool:
        dx = abs(self.x - other.x)
        dy = abs(self.y - other.y)
        dz = abs(self.z - other.z)

        return dx + dy + dz == 1

    def touch_direction(self, other: 'Node'):
        if not self.is_neightbour(other):
            raise NotNeighbourError

        touch_x = abs(self.x - other.x) == 1
        touch_y = abs(self.y - other.y) == 1
        touch_z = abs(self.z - other.z) == 1

        return touch_x, touch_y, touch_z

    def bond_value(self, other: 'Node'):
        return _bond_values.get(frozenset({self.letter, other.letter}), 0)

    def pos_from_direction(self, direction: int, prev: Optional['Node'] = None):
        # Use self.prev if no previous node was specified
        prev = prev if prev else self.prev

        x = prev.x
        y = prev.y
        z = prev.z

        match direction:
            case definitions.UP:
                y += 1
            case definitions.DOWN:
                y -= 1
            case definitions.LEFT:
                x -= 1
            case definitions.RIGHT:
                x += 1
            case definitions.FORWARD:
                z += 1
            case definitions.BACKWARD:
                z -= 1

        return x, y, z

    def change_direction(self, direction: int):
        if not self.direction_from_previous:
            raise Exception("First node in chain!")

        new_x, new_y, new_z = self.pos_from_direction(direction)
        dx = new_x - self.x
        dy = new_y - self.y
        dz = new_z - self.z

        self.x = new_x
        self.y = new_y
        self.z = new_z

        if self.next:
            self.next.cascade_position(dx, dy, dz)

    def cascade_position(self, dx, dy, dz):
        self.x += dx
        self.y += dy
        self.z += dz

        if self.next is not None:
            self.next.cascade_position(dx, dy, dz)
