from experiments_helper import create_experiment_folders, generate_random_sequence
from protein_folding.protein import Protein
from protein_folding.algorithms import DepthFirst
from random import shuffle
import cProfile


def main():
    create_experiment_folders()
    sequence = generate_random_sequence(12)
    score = 0
    dim = 2
    protein = Protein(sequence)

    algorithm = DepthFirst(protein, dimensions=dim, max_iterations=5000, debug=True)
    score = algorithm.run()
    print(f"Sequence: {sequence}")
    print(f"Score: {score}")

    if dim == 2:
        protein.plot(f'./output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{dim}.png')
    elif dim == 3:
        protein.plot_3d(f'./output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{dim}.png')


if __name__ == '__main__':
    # profiler = cProfile.Profile()
    # profiler.runcall(main)
    # profiler.dump_stats('output/protein.prof')

    main()
