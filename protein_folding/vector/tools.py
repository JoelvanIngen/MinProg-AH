from . import Vec3D


def get_min_max(vectors: list[Vec3D]) -> tuple[Vec3D, Vec3D]:
    """
    Creates and returns a tuple consisting of two Vec3D vectors,
        where the first represents the minimum and the second
        the maximum coordinates
    """

    x = [v.x for v in vectors]
    y = [v.y for v in vectors]
    z = [v.z for v in vectors]

    return (Vec3D(min(x), min(y), min(z)),
            Vec3D(max(x), max(y), max(z)))
