from . import Algorithm

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from protein_folding.protein import Protein


class SimulatedAnnealing(Algorithm):
	def __init__(self, protein: 'Protein', dimensions, **kwargs):
		super().__init__(protein, dimensions, **kwargs)

	
