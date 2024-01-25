from tqdm import tqdm

from . import Algorithm
from .heuristics import *
from .pruning import *
from protein_folding.fast_protein import fast_compute_bond_score, fast_validate_protein

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from protein_folding.protein import Protein

direction_dict = {-1: "Left", 1: "Right", 2: "Up", -2: "Down", 3: "Forward", -3: "Backward"}


class DepthFirst(Algorithm):
    """
    A greedy algorithm that applies a direction for each individual
        node, starting from the first node after the root node. It determines
        for each node what the available directions are and which directions
        are already occupied. Then it chooses the direction which yields the
        highest immediate score change. If there are no directions left i.e., it worked
        itself into a corner, it will restart from scratch.
    """

    def __init__(self, protein: 'Protein', dimensions: int, max_iterations: int = 5000, pruning=True, **kwargs):
        super().__init__(protein, dimensions, **kwargs)

        self._iteration = 0
        self.max_iterations = max_iterations
        self.lowest_callback = len(protein)

        self.best_score: int = 1
        self.best_order: list[int] = []

        self.end_nodes_evaluated = 0
        self.amount_of_best_found = 0

        self.dimensions = dimensions

        self.pruning = pruning
        self.budget = 100
        self.cost_per_iteration = max_iterations / dimensions ** (len(self.protein))

    def saving_by_cutting_branch(self, depth: int):
        """
        Calculates the amount of nodes that have (theoretically) been pruned,
            which can then be added back to the budget multiplied with the
            cost per iteration.
        """
        return (2 * self.dimensions - 1) ** (len(self.protein) - depth - 1)

    def _increment_iteration(self):
        self.budget -= 1
        self._iteration += 1

        if self.pbar:
            self.pbar.update(1)

        if self.keep_score_history:
            self.best_score_history.append(self.best_score)

        if self.keep_order_history:
            if fast_validate_protein(self.protein.order[1:]):
                self.order_history.append(self.protein.order[1:])

    def prune(self, directions, depth):
        """
        Prunes the worst branches until either we have enough budget, or until
            there is only one branch left

        Only prune if the following conditions are true:
        - Pruning is enabled
        - We have negative budget
        - There are 2 or more branches remaining
        - We are not at the last iteration (pruning won't save any calculating)
        """
        # print(f"Budget: {self.budget}")
        while self.pruning and self.budget < 0 and len(directions) > 1 and len(self.protein) - depth > 1:
            # Prune worst branch
            directions.pop(-1)
            self.budget += int(self.saving_by_cutting_branch(depth) * self.cost_per_iteration) + 1  # Round up
            # print(f"Nodes pruned: {self.saving_by_cutting_branch(depth)}, new budget: {self.budget} (depth {depth})")

        return directions

    def _next_fold(self, depth: int):
        # Check if we've reached the end
        if depth >= len(self.protein):
            self.end_nodes_evaluated += 1

            bond_score = fast_compute_bond_score(self.protein.sequence, self.protein.order[1:])
            if bond_score < self.best_score:
                self.best_score = bond_score
                self.best_order = self.protein.order
                self.amount_of_best_found = 1
            elif bond_score == self.best_score:
                self.amount_of_best_found += 1

            return

        # # Check if we can prune position
        # pruner = Neighbours(self.protein)
        # if pruner.run() and depth > 1:
        #     # print("SNAP")
        #     return

        pruner = Score(self.protein)
        if pruner.run(best_score=self.best_score, depth=depth):
            # print("SNAP")
            return

        self._increment_iteration()

        if self._iteration > self.max_iterations:
            return

        free_directions = self.protein.nodes[depth].get_free_directions(self.directions)
        if not free_directions:
            return

        direction_scores, free_directions_sorted = self._process_heuristics(
            depth, free_directions)

        free_directions_sorted = self.prune(free_directions_sorted, depth)

        for i, direction in enumerate(free_directions_sorted):
            self.protein.preserve()
            self.protein.nodes[depth].change_direction(direction)
            self._next_fold(depth + 1)
            self.protein.revert()

            # Update lowest node to see how much of the configuration has been tried
            if depth < self.lowest_callback:
                self.lowest_callback = depth

        if self._debug and self._iteration % 100 == 0:
            print(f'Iteration: {self._iteration}/{self.max_iterations}, Depth/lowest: {depth}/{self.lowest_callback},'
                  f' Best score: {self.best_score} ({self.amount_of_best_found} found), End nodes reached: '
                  f'{self.end_nodes_evaluated}')

    def run(self) -> float:
        if self._show_progress:
            self.pbar = tqdm(range(self.max_iterations))

        # Start at first node after root node
        self._next_fold(depth=1)

        if self._debug:
            print(self.best_order)

        self.protein.set_order(self.best_order[1:])
        return self.best_score
