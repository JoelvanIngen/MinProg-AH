import numpy as np
import matplotlib.pyplot as plt

from numpy.typing import NDArray


class Protein:

	def __init__(self, sequence: str, order: NDArray = None) -> None:
		# sequence as ints: 1 = H, 2 = P
		self.sequence = sequence
		#self.seq_as_int = [1 if acid == 'H' else 2 for acid in sequence]
		# initial order: straight line
		if order is not None:
			self.order = order
		else:
			self.order = np.ones(len(sequence) - 1)

	def plot(self) -> None:
		fig = plt.figure()
		plt.xlim(-len(self.sequence), len(self.sequence))
		plt.ylim(-len(self.sequence), len(self.sequence))
		x, y = 0, 0
		for i, acid in enumerate(self.sequence):
			x_new = np.count_nonzero(self.order[:i] == 1) - np.count_nonzero(self.order[:i] == -1)
			y_new = np.count_nonzero(self.order[:i] == 2) - np.count_nonzero(self.order[:i] == -2)
			plt.text(x_new, y_new, acid, size='10')
			plt.plot([x, x_new], [y, y_new], '--', color='black', linewidth=1)
			x, y = x_new, y_new
		plt.savefig("test.png")


if __name__ == "__main__":
	seq = "HHPHHHPH"
	order = np.array([1,2,1,-2,-2,-1,-2])
	protein = Protein(seq, order)
	protein.plot()
