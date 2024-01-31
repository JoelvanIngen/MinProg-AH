import random
from tqdm import tqdm
import sys

from protein_folding.fast_protein import fast_validate_protein, fast_compute_bond_score
from protein_folding.protein import Protein, _get_posvecs_from_order
from protein_folding.vector import Vec3D, get_min_max

from protein_folding.algorithms import Algorithm
from protein_folding.algorithms.heuristics.minimise_dimensions import MinimiseDimensions


class SimulatedAnnealingMinDim(Algorithm):
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
        self.n_permutations = 11000
        self.orders = list()
        self.orders.append(([1] * (len(protein.sequence) - 1)))

    def get_size_score(self, order):
        coords = _get_posvecs_from_order(order)
        left_upper, right_lower = get_min_max(coords)
        box = right_lower - left_upper
        return box.len_sq()

    def get_permutated_directions(self, node_idx: int):
        """
        A function to obtain a set of directions that differs from the current
        protein ordering starting from a direction change at node [node_idx].
        The returned directions are an ordering of the same sequence as
        self.protein, but not necessarily legal.

        Pre:
            - node_idx is an index between 1 and len(sequence) - 1
        Post:
            - returns a list of directions that can be fed to protein.set_order
        """
        dirs_total = [None] + self.protein.get_order()
        node = self.protein.nodes[node_idx]
        free_directions = node.get_free_directions(self.directions)

        if free_directions:
            size_scores = list()
            for direction in free_directions:
                dirs_total[node_idx] = direction 
                size_scores.append(1/self.get_size_score(dirs_total))
                # todo: check validity before this, else smaller but invalid happens often
            
            direction = random.choices(free_directions, weights=size_scores, k=1)[0]
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
        best_order = []
        score = 0
        size_score = 9999
        threshold = 1
        for _ in range(self.n_permutations): #tqdm(range(self.n_permutations)):
            # get random node idx to permutate from and new directions
            node_idx = random.randint(1, len(self.protein.sequence) - 1)
            dirs_total = self.get_permutated_directions(node_idx)

            # Prevent computing score if order is not valid
            if fast_validate_protein(dirs_total):

                comparison_score = fast_compute_bond_score(self.protein.sequence, dirs_total)
                comparison_size_score = self.get_size_score(dirs_total)
                
                decision_float = random.random()
                if (
                    comparison_score < score or
                    (comparison_size_score < size_score and comparison_score == score) or
                    threshold > decision_float
                ):
                    if self.orders[-1] != dirs_total:
                        self.orders.append(dirs_total)
                    self.protein.set_order(dirs_total)
                    score = comparison_score
                    size_score = comparison_size_score

                if comparison_score <= score:
                    best_order = dirs_total

            threshold *= self.decrease

        self.protein.set_order(best_order)
        score = self.protein.get_bond_score()
            
        return score
