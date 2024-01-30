from experiments_helper import import_fix
from tqdm import tqdm
from protein_folding.protein import Protein
from protein_folding.algorithms import *

N_ITERATIONS = 10000
N_DIMENSIONS = 2
SEQUENCE = "HHPCHHPCCPCPPHHHHPPHCHPHPHCHPP"
ALGORITHM = IterativeRandom  # PureRandom or IterativeRandom


def main():
    assert 2 <= N_DIMENSIONS <= 3, f"Number of dimensions must be 2 or 3, but is set to {N_DIMENSIONS}"

    best_score = 1

    for _ in tqdm(range(N_ITERATIONS)):
        p = Protein(SEQUENCE)
        a = ALGORITHM(p, N_DIMENSIONS)
        score = a.run()

        if score < best_score:
            best_score = score

    print(f"Best score: {best_score}")


if __name__ == '__main__':
    main()
