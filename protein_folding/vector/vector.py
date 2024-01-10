from math import sqrt


class Vec3D:
    def __init__(self, x: int, y: int, z: int):
        self.x = x
        self.y = y
        self.z = z

    def __abs__(self) -> 'Vec3D':
        return Vec3D(
            x=abs(self.x),
            y=abs(self.y),
            z=abs(self.z)
        )

    def __add__(self, other: 'Vec3D') -> 'Vec3D':
        return Vec3D(
            x=self.x + other.x,
            y=self.y + other.y,
            z=self.z + other.z
        )

    def __iadd__(self, other: 'Vec3D') -> None:
        self.x += other.x
        self.y += other.y
        self.z += other.z

    def __iter__(self) -> iter:
        return iter((self.x, self.y, self.z))

    def __len__(self) -> float:
        return sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def __sub__(self, other: 'Vec3D') -> 'Vec3D':
        return Vec3D(
            x=self.x - other.x,
            y=self.y - other.y,
            z=self.z - other.z
        )

    @classmethod
    def abs_diff(cls, vec1: 'Vec3D', vec2: 'Vec3D') -> 'Vec3D':
        return abs(vec1 - vec2)

    def area(self) -> int:
        return self.x * self.y

    def eq_components(self, eq_val: int) -> tuple[bool, bool, bool]:
        return self.x == eq_val, self.y == eq_val, self.z == eq_val

    def sum_components(self) -> int:
        return self.x + self.y + self.z

    def len_sq(self) -> int:
        return self.x ** 2 + self.y ** 2 + self.z ** 2

    def volume(self) -> int:
        return self.x * self.y * self.z
