"""
A stack class to save a protein state and load it using as little memory and prevent as many object creations as
possible

Authors: Joël
"""

class Stack:
    def __init__(self):
        self._direction_history: list[list[int | None]] = []
        self._position_dict_history: list[list[dict]] = []

    def push(self, order: list[int | None], pos_dict: list[dict]) -> None:
        self._direction_history.append(order[:])
        self._position_dict_history.append(pos_dict)

    def pull(self) -> tuple[list[int | None], list[dict]]:
        return (
            self._direction_history.pop(),
            self._position_dict_history.pop(),
        )
