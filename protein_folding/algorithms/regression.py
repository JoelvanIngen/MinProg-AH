import random
import pdb
from typing import TYPE_CHECKING
from protein_folding.protein import Protein

from . import Algorithm


class Regression(Algorithm):
    """
    A simple regression algorithm that compares the current protein ordering
    against a randomly permutated but legal state and sets the current state to
    said permutation if the bond score of that state is lower.
    """
    
    def __init__(self, protein: 'Protein', dimensions, **kwargs):
        super().__init__(protein, dimensions, **kwargs)
        # number of random mutations to allow
        self.n_permutations = 200

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
        Runs a simple regression algorithm for n_permutations iterations.
        
        Post:
            - self.protein is ordered in the way that the algorithm found to
            maximise the bond score.
        """
        score = 0
        for _ in range(self.n_permutations):
            # get random node idx to permutate from and new directions
            node_idx = random.randint(1, len(self.protein.sequence) - 1)
            dirs_total = self.get_permutated_directions(node_idx)

            # generate comparison protein to validate score TODO: ugly solution
            comparison_protein = Protein(self.protein.sequence)
            comparison_protein.set_order(dirs_total)
            comparison_score = comparison_protein.get_bond_score()

            # update score and self.protein
            if (comparison_protein.has_valid_order() and
                comparison_score <= self.protein.get_bond_score()):
                self.protein.set_order(dirs_total)
                score = comparison_score

        return score
