import random
import pdb
from typing import TYPE_CHECKING
from protein_folding.fast_protein import fast_validate_protein, fast_compute_bond_score
from protein_folding.protein import Protein

from . import Algorithm


class SimulatedAnnealing(Algorithm):
    """
    A simulated annealing algorithm that compares the current protein ordering
    against a randomly permutated but legal state and sets the current state to
    said permutation if the bond score of that state is lower, and if a
    randomly chosen float is above a perpetually lowering threshold, inspired
    by the way metals anneal.
    """
    
    def __init__(self, protein: 'Protein', dimensions, **kwargs):
        super().__init__(protein, dimensions, **kwargs)
        # decrease of threshold value per iteration
        self.decrease = .999
        # number of random mutations to allow
        # TODO: decide how to determine when to stop algorithm
        self.n_permutations = 500000

    def get_permutated_directions(self, node_idx: int):
        """
        TODO: should be implemented within protein
        A function to obtain a set of directions that differs from the current
        protein ordering starting from a direction change at node [node_idx].
        The returned directions are an ordering of the same sequence as
        self.protein, but not necessarily legal.

        Pre:
            - node_idx is an index between 1 and len(sequence) - 1
        Post:
            - returns a list of directions that can be fed to protein.set_order
        """
        dirs_total = [
            node.direction_from_previous for node in self.protein.nodes
            ]
        node = self.protein.nodes[node_idx]
        free_directions = node.get_free_directions(self.protein.node_positions, self.directions)

        if free_directions:
            direction = random.choice(free_directions)
            dirs_total[node_idx] = direction 
        
        dirs_total = dirs_total[1:]

        return dirs_total
                
    def run(self) -> float:
        """
        Runs a simulated annealing algorithm for n_permutations iterations.
        
        Post:
            - self.protein is ordered in the way that the algorithm found to
            maximise the bond score.
        """
        score = 0
        threshold = 1
        for _ in range(self.n_permutations):
            # get random node idx to permutate from and new directions
            node_idx = random.randint(1, len(self.protein.sequence) - 1)
            dirs_total = self.get_permutated_directions(node_idx)

            # Prevent computing score if order is not valid
            if fast_validate_protein(dirs_total):

                comparison_score = fast_compute_bond_score(self.protein.sequence, dirs_total)
                decision_float = random.random()
                if comparison_score <= score or threshold > decision_float:
                    self.protein.set_order(dirs_total)
                    score = comparison_score

            threshold *= self.decrease
            
        return score
