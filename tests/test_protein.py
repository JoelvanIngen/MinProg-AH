import pytest

import numpy as np
from protein_folding.protein import Protein, InvalidSequenceError


def test_create_valid_protein():
    seq = "HHPHHHPH"
    order = np.array([1, 2, 1, -2, -2, -1, -2])
    _ = Protein(seq, order)


def test_create_invalid_protein():
    with pytest.raises(InvalidSequenceError):
        seq = "HHPHHAPH"
        order = np.array([1, 2, 1, -2, -2, -1, -2])
        _ = Protein(seq, order)
