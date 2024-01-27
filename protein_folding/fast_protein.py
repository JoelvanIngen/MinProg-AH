"""
Contains functions for performance-crucial functionality that perform much better than their corresponding methods from
Protein class.

Author: Joël
"""

import numpy as np

from typing import Sequence


# Map directions to coordinate deltas by getting direction integers and add
# 3 to each index so the list starts at index 0
# Should have fast lookup times thanks to this method
_direction_to_delta = [
	[0, 0, -1],  # Backwards
	[0, -1, 0],  # Down
	[1, 0, 0],   # Right
	[],          # Stop (not implemented)
	[-1, 0, 0],  # Left
	[0, 1, 0],   # Up
	[0, 0, 1],   # Forward
]


def _get_bond_score(letter1: str, letter2: str):
	if letter1 == 'P' or letter2 == 'P':
		return 0

	if letter1 == 'H' and letter2 == 'H':
		return -1
	elif letter1 == 'C' and letter2 == 'C':
		return -5
	else:
		# C-H or H-C, only options left so avoid more comparisons for speed
		return -1


def fast_compute_bond_score(seq: str, order: Sequence[int]) -> int:
	"""
	Computes the bond score for a protein sequence and order. It uses an
		alternative representation that should be much quicker, but harder to
		expand in the future. However, this is not a problem since it will
		only perform a very specific and small-scoped task.
	This alternative representation uses a dictionary to map node positions to
		their letter. During the creation of this dictionary, we count the
		score that is made up by nodes that are next to each other on the
		protein chain. We do this to compensate for over-counting when we count
		the number of nodes that touch each other, since nodes are
		indistinguishable except for their letter. When counting, we only look
		in specific directions to avoid counting twice (Right, Up and Forward).

	pre:
		- seq is a str containing the molecule letters
		- order is a list of ints containing the shape

	post:
		- an int is returned representing the bond value of this order
	"""

	# Keep track of protein score
	score = 0

	# Create dictionary mapping node positions to their letter and add first
	# node
	x = y = z = 0

	nodes = {(x, y, z): seq[0]}

	seq_idx = 1
	order_idx = 0
	for _ in range(len(seq) - 1):
		new_node_delta = _direction_to_delta[order[order_idx] + 3]
		x += new_node_delta[0]
		y += new_node_delta[1]
		z += new_node_delta[2]

		# Do not save node or check score if it's a 'P', since it will not
		# affect the final score
		if seq[seq_idx] != 'P':
			nodes[(x, y, z)] = seq[seq_idx]

			# Remove score to compensate for over-counting in advance
			score -= _get_bond_score(seq[seq_idx], seq[seq_idx - 1])

		seq_idx += 1
		order_idx += 1

	# Go through all values in dict and add score if there is a node right, up
	# or forward
	for pos, letter in nodes.items():
		# Check right
		try_pos = (pos[0] + 1, pos[1], pos[2])
		neighbour_letter = nodes.get(try_pos)
		if neighbour_letter:
			score += _get_bond_score(letter, neighbour_letter)

		# Check up
		try_pos = (pos[0], pos[1] + 1, pos[2])
		neighbour_letter = nodes.get(try_pos)
		if neighbour_letter:
			score += _get_bond_score(letter, neighbour_letter)

		# Check forward
		try_pos = (pos[0], pos[1], pos[2] + 1)
		neighbour_letter = nodes.get(try_pos)
		if neighbour_letter:
			score += _get_bond_score(letter, neighbour_letter)

	return score


def fast_validate_protein(order: Sequence[int]):
	"""
	Computes whether a specified order would result in a valid protein.
	To do this, it uses an alternative representation for the positions of each
	molecule, where we multiply each coordinate by a factor depending on the
	length of the order list, such that we can add each component multiplied to
	this factor and put it in an array, such that each coordinate that can
	theoretically be reached by this order list has its own space, while
	keeping the length of the list at N^3.
	"""

	y_mult = len(order) + 1
	x_mult = y_mult * y_mult + 1

	taken_positions = np.zeros(y_mult ** 3)

	x = y = z = 0
	taken_positions[x * x_mult + y * y_mult + z] = 1

	for direction in order:
		new_node_delta = _direction_to_delta[direction + 3]
		x += new_node_delta[0]
		y += new_node_delta[1]
		z += new_node_delta[2]

		idx = x * x_mult + y * y_mult + z
		# Check if new position is already taken
		if taken_positions[idx] == 1:
			return False

		taken_positions[idx] = 1

	return True
