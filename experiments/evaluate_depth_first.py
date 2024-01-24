from experiments_helper import create_experiment_folders, generate_random_sequence, generate_realistic_sequence
from protein_folding.protein import Protein
from protein_folding.algorithms import DepthFirst
from random import shuffle
import cProfile


def main():
    create_experiment_folders()
    sequence = generate_realistic_sequence(18)
    dim = 2
    protein = Protein(sequence)

    algorithm = DepthFirst(protein, dimensions=dim, max_iterations=20000, debug=True, keep_score_history=True)
    score = algorithm.run()
    print(f"Sequence: {sequence}")
    print(f"Score: {score}")

    algorithm.plot_score_progress()

    if dim == 2:
        protein.plot(f'./output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{dim}.png')
    elif dim == 3:
        protein.plot_3d(f'./output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{dim}.png')


if __name__ == '__main__':
    # profiler = cProfile.Profile()
    # profiler.runcall(main)
    # profiler.dump_stats('output/protein.prof')

    main()
