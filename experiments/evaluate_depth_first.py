"""
Runs the Depth-First algorithm and can tune its parameters if configured to do so.
Author: Joël

Parameters:
    N_ITERATIONS: Amount of iterations the algorithm is allowed to run at most. It might still run fewer iterations if
        pruning is strict.
    N_DIMENSIONS: Can be set to either 2 or 3, and represents the amount of dimensions the protein can fold in.
    USE_RANDOM_SEQUENCE: Determines whether a random sequence will be generated or the custom sequence is used.
    RANDOM_SEQUENCE_LENGTH: Length of sequence to generate if USE_RANDOM_SEQUENCE is set to True.
    CUSTOM_SEQUENCE: Custom sequence to use if USE_RANDOM_SEQUENCE is set to False.
    PLOT_BEST: Visualises the best protein configuration if set to True.
    ALPHA: Alpha parameter to use for pruning. Affects required score threshold increase per depth. ALPHA = 0.5 is a
        good starting value.
    BETA: Beta parameter to use for pruning. Affects fixed score threshold margin. BETA = 18 is a good starting value.
    FIND_BEST_BETA_VALUE: Will try different parameters for Beta using fewer iterations and use the best for full run.
        Very slow but finds good results.
    USE_PROFILING: Switches profiling on or off.

Finding the best value for beta:
    Use FIND_BEST_BETA_VALUE = True and set ALPHA = 0.5. The script will then loop through different Beta values and
    select the one that yields the best results.
    Even better results can be found by manually setting Beta to the best value that the script found, increasing the
    number of iterations and setting Alpha to 0.45 or even 0.40.
"""

import cProfile
from experiments_helper import create_experiment_folders, generate_realistic_sequence
from protein_folding.protein import Protein
from protein_folding.algorithms import DepthFirst
from protein_folding.algorithms.heuristics import *

N_ITERATIONS = 100000
N_DIMENSIONS = 2

USE_RANDOM_SEQUENCE = False
RANDOM_SEQUENCE_LENGTH = 30
CUSTOM_SEQUENCE = "HHPCHHPCCPCPPHHHHPPHCHPHPHCHPP"
PLOT_BEST = True
ANIMATE_SEARCH = False

ALPHA = 0.50
BETA = 18

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

    if FIND_BEST_BETA_VALUE:
        beta = find_best_prune_parameter(sequence)
    else:
        beta = BETA
    print(f"Using alpha = {ALPHA}, beta = {beta} for score-based pruning")

    protein = Protein(sequence)
    algorithm = DepthFirst(protein, dimensions=N_DIMENSIONS, max_iterations=N_ITERATIONS,
                           prune_alpha=ALPHA, prune_beta=beta, debug=True, keep_score_history=True,
                           keep_order_history=ANIMATE_SEARCH, show_progressbar=False,
                           heuristics=[
                               MinimiseDimensions
                           ])

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

    # Prints best configuration depth vs score
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

        if parameter > 1.5 * -best_score:
            parameter = int(1.5 * -best_score)
        else:
            parameter -= 1

    return best_p


if __name__ == '__main__':
    if USE_PROFILING:
        profiler = cProfile.Profile()
        profiler.runcall(main)
        profiler.dump_stats('output/protein.prof')
    else:
        main()
