import sys
import os

sys.path.append('..')


def create_experiment_folders():
    # Ensures folders are created so experiments don't crash
    if not os.path.exists('./output/'):
        os.makedirs('./output/')


def create_bruteforce_folders():
    # Ensures Brute Force folders are created so experiments don't crash
    if not os.path.exists('./bf_output/'):
        os.makedirs('./bf_output/')
