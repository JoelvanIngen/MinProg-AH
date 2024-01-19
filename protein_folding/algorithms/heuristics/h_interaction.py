
from . import Heuristic
from protein_folding.vector import Vec3D

class HydrophobicInteraction(Heuristic):
    def __init__(self, *args):
        super().__init__(*args)
        hydrophobic_letters = {'H'}
        self.hydrophobic_atoms = [node for node in self.protein.nodes if node.letter in hydrophobic_letters]

    def _calculate_potential(self, atom):
        # Compute the potential for a single hydrophobic atom
        potential = 0
        for other_atom in self.hydrophobic_atoms:
            if other_atom != atom:
                distance = Vec3D.distance(atom.position, other_atom.position)
                potential -= 1 / distance if distance != 0 else float('inf')
        return potential

    def run(self):
        total_potential = 0
        for atom in self.hydrophobic_atoms:
            total_potential += self._calculate_potential(atom)
        return total_potential
