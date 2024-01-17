from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms.hillclimb import HillClimb


# Create a Protein instance with a specific sequence or configuration
protein_instance = Protein(sequence="HHPHHHPPH")

# Initialize the Hill Climb algorithm
hill_climb = HillClimb(protein=protein_instance, dimensions = 2, temperature=1000.0, cooling_rate=0.001)

# Run the algorithm
optimized_protein = hill_climb.run()
