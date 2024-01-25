from experiments_helper import create_experiment_folders, generate_random_sequence, generate_realistic_sequence
from protein_folding.protein import Protein
from protein_folding.algorithms import DepthFirst
from protein_folding.algorithms.heuristics import *
from random import shuffle
import cProfile


def main():
    create_experiment_folders()
    sequence = generate_realistic_sequence(30)
    # sequence = "HHPCHHPCCPCPPHHHHPPHCHPHPHCHPP"
    dim = 2
    protein = Protein(sequence)

    algorithm = DepthFirst(protein, dimensions=dim, max_iterations=20000, prune_multiplier=7, debug=True, keep_score_history=True,
                           keep_order_history=True,
                           heuristics=[
                               # PotentialPlus,
                               MinimiseDimensions
                           ])
    score = algorithm.run()
    print(f"Sequence: {sequence}")
    print(f"Score: {score}")

    if dim == 2:
        protein.plot(f'./output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{dim}.png')
    elif dim == 3:
        protein.plot_3d(f'./output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{dim}.png')

    # algorithm.plot_score_progress()
    protein.animate_2d(algorithm.order_history)

    from protein_folding.fast_protein import fast_compute_bond_score
    for i in range(4, len(protein.order) - 1):
        print(i, fast_compute_bond_score(protein.sequence[:i], protein.order[1:i]))


if __name__ == '__main__':
    # profiler = cProfile.Profile()
    # profiler.runcall(main)
    # profiler.dump_stats('output/protein.prof')

    main()
