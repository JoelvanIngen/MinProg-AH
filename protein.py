import numpy as np
import matplotlib.pyplot as plt


class Protein:

	def __init__(self, sequence: str, order: list = None) -> None:
		# sequence as ints: 1 = H, 2 = P
		self.sequence = sequence
		#self.seq_as_int = [1 if acid == 'H' else 2 for acid in sequence]
		# initial order: straight line
		if order:
			self.order = order
		else:
			self.order = np.ones(len(sequence) - 1)

	def plot_grid(self) -> None:
		# TODO


if __name__ == "__main__":
	seq = "abcdef"
	order = [1,2,1,-2,1]
	protein = Protein(seq, order)
