from . import Heuristic
from protein_folding.protein import Protein
from protein_folding.node import Node
from protein_folding.vector.tools import get_min_max  

class MinimizeDimensions(Heuristic):
    def __init__(self, *args):
        super().__init__(*args)

    def _calc_size_score(self):
        """
        Computes a score based on the physical dimensions of the protein. This
            value can be used as a penalty for a specific configuration.
        """
        dim = get_min_max([node.pos for node in self.nodes])
        box = dim[1] - dim[0]
        return box.len_sq()
