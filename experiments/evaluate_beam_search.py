import cProfile
from experiments_helper import create_experiment_folders, generate_realistic_sequence
from protein_folding.protein import Protein
from protein_folding.algorithms import BeamSearch
from protein_folding.algorithms.heuristics import *

N_ITERATIONS = 100000
N_DIMENSIONS = 2

USE_RANDOM_SEQUENCE = False
RANDOM_SEQUENCE_LENGTH = 30
CUSTOM_SEQUENCE = "HHPHHHPHPHHHPH"
PLOT_BEST = True
ANIMATE_SEARCH = False

FIND_BEST_BETA_VALUE = True
USE_PROFILING = False


def main():
    assert 2 <= N_DIMENSIONS <= 3, f"Number of dimensions must be 2 or 3, but is set to {N_DIMENSIONS}"

    create_experiment_folders()

    if USE_RANDOM_SEQUENCE:
        sequence = generate_realistic_sequence(RANDOM_SEQUENCE_LENGTH)
    else:
        sequence = CUSTOM_SEQUENCE
    print(f"Using sequence {sequence}")

    protein = Protein(sequence)
    algorithm = BeamSearch(protein, dimensions=N_DIMENSIONS, max_iterations=N_ITERATIONS,
                             debug=True, keep_score_history=True,
                             keep_order_history=ANIMATE_SEARCH, show_progressbar=True)

    score = algorithm.run()
    print(f"Sequence: {sequence}")
    print(f"Score: {score}")

    if PLOT_BEST:
        filename = f'./output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{N_DIMENSIONS}.png'
        plotter = protein.plot if N_DIMENSIONS == 2 else protein.plot_3d

        plotter(filename)

    # algorithm.plot_score_progress()

    if ANIMATE_SEARCH:
        protein.animate_2d(algorithm.order_history)


if __name__ == '__main__':
    if USE_PROFILING:
        profiler = cProfile.Profile()
        profiler.runcall(main)
        profiler.dump_stats('output/protein.prof')
    else:
        main()
