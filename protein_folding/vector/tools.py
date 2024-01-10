from . import Vec3D


def get_min_max(vectors: list[Vec3D]) -> Vec3D:
    x = [v.x for v in vectors]
    y = [v.y for v in vectors]
    z = [v.z for v in vectors]

    dx = max(x) - min(x)
    dy = max(y) - min(y)
    dz = max(z) - min(z)

    return Vec3D(dx, dy, dz)
