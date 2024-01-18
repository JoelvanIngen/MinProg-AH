import sys
sys.path.append('../../..')

from protein_folding.algorithms.bruteforce import BruteForce
from sequence_generator import sequence_generator

def main():
	print(sequence_generator(50))


if __name__ == "__main__":
	main()