import sys
from protein_folding.protein import Protein
from protein_folding.algorithms import *
from protein_folding.algorithms.heuristics import *


def get_user_int(*args):
    while True:
        user_string = input(*args)

        try:
            return int(user_string)
        except ValueError:
            print("Please input a number.")


def get_user_float(*args):
    while True:
        user_string = input(*args)

        try:
            return float(user_string)
        except ValueError:
            print("Please input a number.")


def get_user_bool(*args):
    while True:
        user_string = input(*args).lower()

        if user_string in ['y', 'yes']:
            return True
        elif user_string in ['n', 'no']:
            return False

        print("Please input y/n.")


def _all_algorithms():
    return {
        'PureRandom': PureRandom,
        'IterativeRandom': IterativeRandom,
        'Spiral': Spiral,
        'BruteForce': BruteForce,
        'Regression': Regression,
        'SimulatedAnnealing': SimulatedAnnealing,
        'Greedy': Greedy,
        'BeamSearch': BeamSearch,
        'DepthFirst': DepthFirst,
    }


def choose_dimensions() -> int:
    while True:
        n_dim = get_user_int("Choose a number of dimensions (2 / 3): ")

        if 2 <= n_dim <= 3:
            return n_dim
        else:
            print("Please input a value of 2 or 3.")


def choose_algorithm():
    algorithms = _all_algorithms()
    algorithm_order = []

    print("Available algorithms:")
    for i, (algorithm_name, algorithm) in enumerate(algorithms.items()):
        print(f"[{i}] {algorithm_name}")
        algorithm_order.append(algorithm)

    while True:
        algorithm_number = get_user_int(f"Choose an algorithm 0-{len(algorithms) - 1}: ")
        if 0 <= algorithm_number <= len(algorithms) - 1:
            return algorithm_order[algorithm_number]


def get_sequence() -> str:
    match len(sys.argv):
        case 1:
            sequence = input("Enter sequence: ")
        case 2:
            sequence = sys.argv[2]
        case _:
            print("Usage: python main.py [sequence] or python main.py and input sequence manually.")
            sys.exit(1)

    return sequence


def set_parameters(algorithm):
    settable_params = algorithm.user_parameters

    for param, default in settable_params:
        change_param = get_user_bool(f"Change parameter {param} (default {default})? (y/n): ")
        if not change_param:
            continue

        value = get_user_float(f"Set parameter {param} to: ")
        setattr(algorithm, param, value)


def set_bool_parameters(algorithm):
    setattr(algorithm, '_debug', get_user_bool("Use debugging? (y/n): "))
    setattr(algorithm, '_visualise', get_user_bool("Keep state history for animating? (y/n): "))
    setattr(algorithm, '_show_progress', get_user_bool("Use progressbar? (y/n): "))


def main():
    sequence = get_sequence()
    protein = Protein(sequence)

    chosen_algorithm = choose_algorithm()
    n_dims = choose_dimensions()

    algorithm = chosen_algorithm(protein, dimensions=n_dims)

    set_parameters(algorithm)

    set_bool_parameters(algorithm)

    score = algorithm.run()
    print(f"\n\nFinal score: {score}")

    # Visualisation??
    protein.plot()


if __name__ == '__main__':
    main()
