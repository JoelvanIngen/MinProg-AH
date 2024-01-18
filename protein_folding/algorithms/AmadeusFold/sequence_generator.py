import random

def sequence_generator(length: int, fractions: list = [.5, .25, .25]) -> list:
	"""
	generate a sequence of len length of valid molecules (letters). Fractions
	denotes the fraction of the total letters that will be (not exactly)
	composed of that letter. Ex: for fractions = [.5, .25, .25], half of the
	letters in the sequence will be H, a quarter will be C and a quarter will
	be P.
	
	pre:
		- length is an int > 0
		- fractions is a list of len = 3 containing floats, the sum of which
		is 1.0
	post:
		- returns a list of length = length made up of H,  C, P
	"""
	sequence = ['H'] * int(length * .5)
	sequence += ['C'] * int(length * .25)
	sequence += ['P'] * (length - int(length * .5) - int(length * .25))
	random.shuffle(sequence)
	return sequence