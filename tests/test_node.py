import pytest

from protein_folding.protein import Protein
from protein_folding.definitions import *


def test_available_directions_all_available():
    seq = "HHPHHHPH"
    p = Protein(seq)

    n = p.nodes[1]

    dirs = [LEFT, RIGHT, UP, DOWN]
    assert n.get_free_directions(p.node_positions, dirs) == [LEFT, RIGHT, UP, DOWN]
