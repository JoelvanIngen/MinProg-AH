import experiments_helper
import time
from protein_folding.protein import Protein

seq = 'H' * 30
order = [1] * 29
N = 50_000


def main():
    time_start = time.perf_counter()
    run()
    time_stop = time.perf_counter()

    time_elapsed = time_stop - time_start

    print(time_elapsed)


def run():
    for _ in range(N):
        create_and_validate_protein()


def create_and_validate_protein():
    p = Protein(seq)
    p.set_order(order)
    score = p.get_bond_score()
    return score


if __name__ == '__main__':
    main()
