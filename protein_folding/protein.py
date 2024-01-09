import numpy as np
import matplotlib.pyplot as plt
from protein_folding.definitions import *
from protein_folding.node import Node


_protein_letter_mapping = {
	'H': 1,
	'P': 2,
	'C': 3,
}


class InvalidSequenceError(Exception):
	pass


class Protein:

	def __init__(self, sequence: str) -> None:
		"""
		representation of a simple protein consisting of polar or hydrophobic
		amino acids in a chain. 

		pre:
			- sequence should be string containing only letter H and P
		post:
			- protein is initialised in straight line ordering
		"""
		# Ensure sequence only consists of valid letters
		if not all(c in _protein_letter_mapping.keys() for c in sequence):
			raise InvalidSequenceError
		self.sequence = sequence
		self.encoded_sequence = np.asarray([_protein_letter_mapping[letter] for letter in sequence], dtype=np.int8)

		# Create list of all the nodes
		self.nodes = [Node(self.sequence[0], 0, 0, 0, direction=None)]

		for c in self.sequence[1:]:
			# Initalise new node in a straight line
			self.nodes.append(Node.from_previous(c, RIGHT, self.nodes[-1]))

		# initial order: straight line to right
		self.set_order(list(np.ones(len(sequence) - 1)))

	def set_order(self, order: list) -> None:
		"""
		sets order (shape) of protein to provided order, and recalculates
		coordinates of all amino acids under new order. Orders are formatted
		as a list of 1, -1, 2, -2 in some order, corresponding to left, right,
		up and down respectively, ie: [1,1,-2,-1,2] (len = len(seq) - 1)

		pre:
			- order is a list formatted as described
		post:
			- attributes for acid coordinates and order are updated
		"""
		# TODO: determine requirements of valid order
		self.order = order 
		# list of tuples containing (x, y)-coordinates for all present acids
		self.acid_coords = list()
		for i in range(len(self.sequence)):
			x = self.order[:i].count(LEFT) - self.order[:i].count(RIGHT)
			y = self.order[:i].count(UP) - self.order[:i].count(DOWN)
			self.acid_coords.append((x, y))
	
	def get_order_quality(self) -> float:
		"""
		quantises the quality of the current fold, or order. Depends on bonds
		and compactness (?)

		pre:
		post:
			- a float is returned representing the fold quality
		"""
		# TODO: quantify compactness (n needlessly empty blocks?)
		n_H_bonds = 0
		for i, crds_1 in enumerate(self.acid_coords):
			for j, crds_2 in enumerate(
					self.acid_coords[:i - 1] + self.acid_coords[i + 2:]):
				if (
					# "nearby": 1 space away
					((crds_1[0] == crds_2[0] and 
						crds_1[1] - crds_2[1] == 1) or
					(crds_1[1] == crds_2[1] and 
						crds_1[0] - crds_2[0]) == 1) and
					# ignore next and previous in chain
					self.sequence[i] == 'H' and
					self.sequence[self.acid_coords.index(crds_2)] == 'H'
				):
					n_H_bonds += 1

		return -n_H_bonds # TODO: compactness

	def plot(self, filename="./unnamed_protein.png") -> None:
		"""
		save current protein configuration as image.
		pre:
			- filename is a string referring to a valid storage location and
			ending in a valid image format
		post:
			- order is saved as image under filename
		"""
		fig = plt.figure()
		x, y = 0, 0
		for i, acid in enumerate(self.sequence):
			x_new, y_new = self.acid_coords[i]
			# place character
			plt.text(x_new, y_new, acid, size='10')
			# place line
			plt.plot([x, x_new], [y, y_new], '--', color='black', linewidth=1)
			x, y = x_new, y_new

		# plot settings
		plt.xlim(-len(self.sequence), len(self.sequence))
		plt.ylim(-len(self.sequence), len(self.sequence))
		plt.axis('off')
		# save image
		plt.savefig(filename)


if __name__ == "__main__":
	seq = "PHPPHPHH"
	order = [2,-1,2,1,1,-2,-2]
	protein = Protein(seq)
	protein.set_order(order)
	print(f"Order quality: {protein.get_order_quality()}") # should be -2.0
	protein.plot()
