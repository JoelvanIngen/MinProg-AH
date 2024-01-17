from protein_folding.protein import Protein
from protein_folding.definitions import *


def test_protein_quality_square():
    seq = "CHHC"
    protein = Protein(seq)
    order = [LEFT, UP, RIGHT]
    protein.set_order(order)
    assert protein.get_bond_score() == -5
