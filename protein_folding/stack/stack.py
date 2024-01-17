class Stack:
    def __init__(self):
        self._order_history: list[list[int]] = []
        self._ghost_history: list[list[bool]] = []

    def push(self, order: list[int], ghost: list[bool]):
        self._order_history.append(order)
        self._ghost_history.append(ghost)
