import experiments_helper
import time
import cProfile
from protein_folding.protein import Protein
from protein_folding.fast_protein import fast_validate_protein

seq = 'H' * 30
order = [1] * 29
N = 5_000


def main():
    # Actual Protein class method
    # profiler1 = cProfile.Profile()
    # profiler1.enable()
    time_start = time.perf_counter()
    res1 = run_protein()
    time_stop = time.perf_counter()
    # profiler1.disable()
    # profiler1.dump_stats('output/protein.prof')

    time_elapsed = time_stop - time_start

    print(f"Time elapsed using Protein class: {round(time_elapsed, 5)}s "
          f"({round(time_elapsed / N * 10_000, 5)}s / 10.000 iterations)")

    # Fast protein method
    # profiler2 = cProfile.Profile()
    # profiler2.enable()
    time_start = time.perf_counter()
    res2 = run_fast()
    time_stop = time.perf_counter()
    # profiler2.disable()
    # profiler2.dump_stats('output/fast_protein.prof')

    time_elapsed = time_stop - time_start

    print(f"Time elapsed using fast function: {round(time_elapsed, 5)}s "
          f"({round(time_elapsed / N * 10_000, 5)}s / 10.000 iterations)")

    if res1 == res2:
        print(f"Sanity check: both methods yielded same results: {res1}")
    else:
        print(f"!!! Methods yielded different results: {res1}, {res2}")


def run_protein():
    result = 0
    for _ in range(N):
        result = create_and_validate_protein()
    return result


def run_fast():
    result = 0
    for _ in range(N):
        result = create_and_validate_fast()
    return result


def create_and_validate_protein():
    p = Protein(seq)
    p.set_order(order)
    score = p.has_valid_order()
    return score


def create_and_validate_fast():
    score = fast_validate_protein(order)
    return score


if __name__ == '__main__':
    main()
