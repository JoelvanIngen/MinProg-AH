import random
import pdb
from tqdm import tqdm
from typing import TYPE_CHECKING
from protein_folding.fast_protein import fast_validate_protein, fast_compute_bond_score
from protein_folding.protein import Protein

from . import Algorithm


class QTable:
    def __init__(self):
        self.table = {}

    def get(self, state, action):
        return self.table.get((state, action), 0)

    def set(self, state, action, value):
        self.table[(state, action)] = value

    def update(self, state, action, reward, learning_rate, discount_factor, next_state, possible_actions):
        old_value = self.get(state, action)
        future_rewards = [self.get(next_state, next_action) for next_action in possible_actions]
        max_future_reward = max(future_rewards, default=0)
        new_value = old_value + learning_rate * (reward + discount_factor * max_future_reward - old_value)
        self.set(state, action, new_value)
