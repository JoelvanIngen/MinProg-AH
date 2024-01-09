_bond_values = {
    {'H', 'H'}: 1,
    {'C', 'C'}: 5,
    {'C', 'H'}: 1,
}


class NotNeighbourError(Exception):
    pass


class Node:
    def __init__(self, letter, x, y, z):
        self.letter = letter
        self.x: int = x
        self.y: int = y
        self.z: int = z

        self._neighbours: list[Node]

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
        return _bond_values.get({self.letter, other.letter}, 0)
