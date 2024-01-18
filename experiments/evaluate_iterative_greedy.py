from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import IterativeGreedy
from random import shuffle
import cProfile


def main():
    create_experiment_folders()
    sequence = 'H' * 3
    sequence += 'C' * 3
    sequence += 'P' * 3
    l = list(sequence)
    shuffle(l)
    sequence = ''.join(l)
    # sequence = "HHPHHHPHPHHHPH"
    score = 0
    dim = 2
    protein = Protein(sequence)

    algorithm = IterativeGreedy(protein, dimensions=dim, max_iterations=100, debug=True)
    score = algorithm.run()
    print(f"Score: {score}")

    if dim == 2:
        protein.plot(f'./output/evaluate_greedy_protein_len{len(sequence)}_dim{dim}.png')
    elif dim == 3:
        protein.plot_3d(f'./output/evaluate_greedy_protein_len{len(sequence)}_dim{dim}.png')


if __name__ == '__main__':
    profiler = cProfile.Profile()
    profiler.runcall(main)
    profiler.dump_stats('output/protein.prof')