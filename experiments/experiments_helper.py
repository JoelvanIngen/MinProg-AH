import os


def create_experiment_folders():
    # Ensures folders are created so experiments don't crash
    if not os.path.exists('./output/'):
        os.makedirs('./output/')
