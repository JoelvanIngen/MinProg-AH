import pytest

import numpy as np
from protein_folding.protein import Protein, InvalidSequenceError


def test_create_valid_protein():
    seq = "HHPHHHPH"
    _ = Protein(seq)


def test_create_invalid_protein():
    with pytest.raises(InvalidSequenceError):
        seq = "HHPHHAPH"
        _ = Protein(seq)
