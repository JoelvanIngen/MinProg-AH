import pytest

import numpy as np
from protein_folding.protein import Protein, InvalidSequenceError
from protein_folding.definitions import *


def test_create_valid_protein():
    seq = "HHPHHHPH"
    _ = Protein(seq)


def test_create_invalid_protein():
    with pytest.raises(InvalidSequenceError):
        seq = "HHPHHAPH"
        _ = Protein(seq)


def test_protein():
    seq = "HHHH"
    protein = Protein(seq)
    order = [LEFT, LEFT, LEFT]
    protein.set_order(order)
    assert protein.get_order_quality() == -3
