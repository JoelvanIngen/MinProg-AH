from . import Heuristic
from protein_folding.protein import Protein
from protein_folding.node import Node
from protein_folding.vec3d import Vec3D  # Assuming Vec3D is a class for 3D vector operations

class MinimizeDimensions(Heuristic):
    def __init__(self, *args):
        super().__init__(*args)


