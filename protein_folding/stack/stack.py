class Stack:
    def __init__(self):
        self._direction_history: list[list[int]] = []
        self._ghost_history: list[list[bool]] = []
        self._position_dict_history: list[dict] = []

    def push(self, order: list[int], ghost: list[bool], pos_dict: dict) -> None:
        self._direction_history.append(order)
        self._ghost_history.append(ghost)
        self._position_dict_history.append(pos_dict)

    def pull(self) -> tuple[list[int], list[bool], dict]:
        return (
            self._direction_history.pop(),
            self._ghost_history.pop(),
            self._position_dict_history.pop(),
        )
