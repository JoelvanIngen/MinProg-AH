from . import Heuristic
from protein_folding.protein import Protein
from protein_folding.node import Node

class foldamount(Heuristic):
    def __init__(self, *args):
        super().__init__(*args)

    def _find_non_ordered_links(self, node_idx: int, direction: int) -> list[tuple]:
        # Method to find non-ordered links for a given node in a specific direction.
        test_protein = Protein(self.protein.sequence)
        test_protein.set_order(self.protein.get_order())
        
        # Change direction of the node
        test_protein.nodes[node_idx].change_direction(direction)

        non_ordered_links = []

        for i, node in enumerate(test_protein.nodes):
            if i != node_idx and abs(i - node_idx) > 1 and test_protein.nodes[node_idx].is_neighbour(node):
                non_ordered_links.append((node, test_protein.nodes[node_idx]))
        # .direction = previousdirection -> 
        return non_ordered_links     
    
    def _calculate_links_score(self, non_ordered_links: list[tuple]) -> float:
        # For now every ordered link is equal!
        return len(non_ordered_links)




