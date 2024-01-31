import heapq
from tqdm import tqdm

from protein_folding.fast_protein import fast_compute_bond_score, fast_validate_protein
from . import Algorithm

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from protein_folding.protein import Protein


class QueueEmptyError(Exception):
    pass


class StateQueueHeapQ:
    def __init__(self):
        self._queue = []
        heapq.heappush(self._queue, ((0, 0), [1]))

        self.iteration = 1

    def get_highest_score_state(self) -> list[int]:
        """
        Finds and returns state with the highest score so far
        """
        if len(self._queue) == 0:
            raise QueueEmptyError

        _score_iteration, state = heapq.heappop(self._queue)
        print(_score_iteration[0])
        return state

    def push(self, order: list[int], score: int) -> None:
        heapq.heappush(self._queue, ((-score, self.iteration), order))
        self.iteration += 1


class StateQueueStandard:
    def __init__(self):
        self.scores: list[int] = [0]
        self.states: list[list[int]] = [[1]]

    def get_highest_score_state(self) -> list[int]:
        """
        Finds and returns state with the highest score so far
        """
        if len(self.scores) == 0:
            assert len(self.states) == 0
            raise QueueEmptyError

        min_score = min(self.scores)
        idx = self.scores.index(min_score)

        self.scores.pop(idx)
        return self.states.pop(idx)

    def push(self, order: list[int], score: int) -> None:
        self.scores.append(score)
        self.states.append(order)


class BeamSearch(Algorithm):
    def __init__(self, protein: 'Protein', dimensions: int, max_iterations: int = 5000, **kwargs):
        super().__init__(protein, dimensions, **kwargs)

        self._iteration = 0
        self.max_iterations = max_iterations

        self.best_score: int = 1
        self.best_order: list[int] = []

        self.end_nodes_evaluated = 0
        self.amount_of_best_found = 0

        self.queue = StateQueueStandard()

        self.user_parameters = [('max_iterations', self.max_iterations)]

    def _increment_iteration(self):
        if self.pbar:
            self.pbar.update(1)

        # UNUSED
        if self.keep_score_history:
            self.best_score_history.append(self.best_score)

        # UNUSED
        if self.keep_order_history:
            if fast_validate_protein(self.protein.order):
                self.order_history.append(self.protein.order)

    def find_next_orders(self, order) -> list[list[int]]:
        valid_next_orders = []
        for direction in self.directions:
            if direction == -order[-1]:
                continue

            test_order = order + [direction]
            if not fast_validate_protein(test_order):
                continue

            valid_next_orders.append(test_order)

        return valid_next_orders

    def get_best_order(self):
        state = self.queue.get_highest_score_state()
        return state

    def save_orders(self, orders):
        for order in orders:
            self.save_order(order)

    def save_order(self, order):
        score = fast_compute_bond_score(self.protein.sequence[:len(order)], order)
        self.queue.push(order, score)

    def keep_best_scores(self, orders):
        for order in orders:
            score = fast_compute_bond_score(self.protein.sequence[:len(order) + 1], order)
            if score < self.best_score:
                if self.verbose:
                    print(f"Found new best score of {score} with order: {order}")

                self.best_score = score
                self.best_order = order[:]

    def process_order(self, order):
        next_orders = self.find_next_orders(order)
        if not next_orders:
            return

        depth = len(next_orders[0])
        if depth == len(self.protein) - 1:
            return self.keep_best_scores(next_orders)

        self.save_orders(next_orders)

        if self.debug and self._iteration % 1000 == 0:
            print(f"Depth: {depth} | "
                  f"Score so far: {fast_compute_bond_score(self.protein.sequence[:len(next_orders[0]) + 1], next_orders[0])}")

    def run(self) -> int:
        if self.show_progress:
            self.pbar = tqdm(range(self.max_iterations))

        # Start at first node after root node
        for iteration in range(self.max_iterations):
            try:
                order = self.get_best_order()
            except QueueEmptyError:
                break

            self.process_order(order)

            self._increment_iteration()

        if self.verbose:
            print(f"Best order: {self.best_order}")

        self.protein.set_order(self.best_order)
        return self.best_score
