import sys
from .protein import Protein


def get_sequence() -> str:
    match len(sys.argv):
        case 1:
            sequence = input("Enter sequence: ")
        case 2:
            sequence = sys.argv[2]
        case _:
            print("Usage: python main.py [sequence] or python main.py and input sequence manually.")
            sys.exit(1)

    return sequence


def main():
    sequence = get_sequence()
    protein = Protein(sequence)

    # Algorithm stuff??
    # algorithm = Algorithm(protein)
    # algorithm.optimise() / algorithm.run() or any cool name we come up with

    # Visualisation??
    protein.plot()

    # TODO: omg vet cool


if __name__ == '__main__':
    main()
