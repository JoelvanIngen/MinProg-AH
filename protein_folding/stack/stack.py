from copy import deepcopy


class Stack:
    def __init__(self):
        self._direction_history: list[list[int | None]] = []
        self._ghost_history: list[list[bool]] = []
        self._position_dict_history: list[dict] = []

    def push(self, order: list[int | None], ghost: list[bool], pos_dict: dict) -> None:
        self._direction_history.append(deepcopy(order))
        self._ghost_history.append(deepcopy(ghost))
        self._position_dict_history.append(deepcopy(pos_dict))

    def pull(self) -> tuple[list[int | None], list[bool], dict]:
        return (
            self._direction_history.pop(),
            self._ghost_history.pop(),
            self._position_dict_history.pop(),
        )
