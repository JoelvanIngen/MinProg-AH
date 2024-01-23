# Ignore Pycharm unused import warnings and keep import_fix!
import import_fix
import random
import os
from protein_folding.algorithms.heuristics import *


def create_experiment_folders():
    # Ensures folders are created so experiments don't crash
    if not os.path.exists('./output/'):
        os.makedirs('./output/')


def create_bruteforce_folders():
    # Ensures Brute Force folders are created so experiments don't crash
    if not os.path.exists('output/bf_output/'):
        os.makedirs('output/bf_output/')


def generate_random_sequence(length):
    return ''.join(random.choices(['H', 'P', 'C'], k=length))


def get_available_heuristics():
    return (
        FoldAmount,
        MinimiseDimensions,
        Potential,
    )
