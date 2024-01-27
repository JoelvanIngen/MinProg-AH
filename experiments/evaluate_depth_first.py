import cProfile
from experiments_helper import create_experiment_folders, generate_random_sequence, generate_realistic_sequence
from protein_folding.protein import Protein
from protein_folding.algorithms import DepthFirst
from protein_folding.algorithms.heuristics import *

N_ITERATIONS = 1000000
N_DIMENSIONS = 2

USE_RANDOM_SEQUENCE = False
RANDOM_SEQUENCE_LENGTH = 30
CUSTOM_SEQUENCE = "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH"

ALPHA = 0.5
BETA = 10

FIND_BEST_PRUNE_PARAMETERS = False
USE_PROFILING = False


def main():
    create_experiment_folders()

    if USE_RANDOM_SEQUENCE:
        sequence = generate_realistic_sequence(RANDOM_SEQUENCE_LENGTH)
    else:
        sequence = CUSTOM_SEQUENCE
    print(f"Using sequence {sequence}")

    if FIND_BEST_PRUNE_PARAMETERS:
        beta = find_best_prune_parameter(sequence)
    else:
        beta = BETA
    print(f"Using beta = {beta} for score-based pruning")

    protein = Protein(sequence)
    algorithm = DepthFirst(protein, dimensions=N_DIMENSIONS, max_iterations=N_ITERATIONS,
                           prune_alpha=ALPHA, prune_beta=beta, debug=True, keep_score_history=True,
                           keep_order_history=True,
                           heuristics=[
                               # PotentialPlus,
                               MinimiseDimensions
                           ])

    score = algorithm.run()
    print(f"Sequence: {sequence}")
    print(f"Score: {score}")

    if N_DIMENSIONS == 2:
        protein.plot(f'./output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{N_DIMENSIONS}.png')
    elif N_DIMENSIONS == 3:
        protein.plot_3d(f'./output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{N_DIMENSIONS}.png')

    # algorithm.plot_score_progress()
    protein.animate_2d(algorithm.order_history)

    from protein_folding.fast_protein import fast_compute_bond_score
    for i in range(4, len(protein.order) - 1):
        print(i, fast_compute_bond_score(protein.sequence[:i], protein.order[1:i]))


def find_best_prune_parameter(sequence: str):
    run_iterations = int(N_ITERATIONS/2)
    parameter = int(len(sequence) * (5/6))
    best_score = 1
    best_p = 15

    while parameter > 2.5:
        print(f"Trying prune_parameter = {parameter}")
        protein = Protein(sequence)
        algorithm = DepthFirst(protein, dimensions=N_DIMENSIONS, max_iterations=run_iterations,
                               prune_alpha=ALPHA, prune_beta=parameter, keep_score_history=True,
                               keep_order_history=True, show_progressbar=False,
                               heuristics=[
                                   MinimiseDimensions
                               ])

        run_score = algorithm.run()
        print(f"Ran algorithm, found score of {run_score} with parameter = {parameter}. "
              f"Ran {algorithm._iteration}/{run_iterations} iterations.")

        if run_score < best_score:
            best_score = run_score
            best_p = parameter

        if algorithm._iteration < run_iterations / 2:
            # More pruning will only result in even less searches
            break

        parameter -= 1

    return best_p


if __name__ == '__main__':
    if USE_PROFILING:
        profiler = cProfile.Profile()
        profiler.runcall(main)
        profiler.dump_stats('output/protein.prof')
    else:
        main()
