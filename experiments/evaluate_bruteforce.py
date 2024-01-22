import cProfile
import pstats

from experiments_helper import create_bruteforce_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import BruteForce
import time


def run():
    # Create necessary folders
    create_bruteforce_folders()

    seq = "HHPCHHPCCPCCCC"
    max_iterations = 0
    dim = 2
    protein = Protein(seq)
    algorithm = BruteForce(protein, dimensions=dim, max_iterations=max_iterations, verbose=True)

    results = algorithm.run()
    print(results)
    print(f"Lowest score: {min(results.values())}")

    for i, order in enumerate(results.keys()):
        protein = Protein(seq)
        protein.set_order(order)
        score = results[order]
        print(f"Protein {i + 1}: {score}")
        protein.plot(filename=f'./output/bf_output/{algorithm.get_name()}_#{i + 1}_score{score}.png')


if __name__ == '__main__':
    # profiler = cProfile.Profile()
    # profiler.runcall(run)
    # stats = pstats.Stats(profiler)
    # stats.strip_dirs().sort_stats('cumulative').print_stats(10)
    # profiler.dump_stats('output/protein.prof')
    start_time = time.perf_counter()
    run()
    print(time.perf_counter() - start_time)
