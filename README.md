# MinProg-AH: Protein Po(w)der
Jarec Schouten, Joël van Ingen & Wolf Gautier

## Experiments (milestone)
To run timed tests on simulated annealing algorithm: run experiments/iterate_sa.py, adjust sequence and timing as necessary.

## Algorithms
**Brute Force**: Tries all possible combinations. Given enough time, this algorithm will find the best shape, but is impractical to run for sufficiently large proteins.

**PureRandom**: Generates a random list of directions, ensures it's valid and computes the resulting score.

**IterativeRandom**: Loops through the protein, and finds the available directions for each node separately, trying to prevent collisions.

**Spiral**: Arranges the protein in a spiral shape.

**Greedy**: Loops through the protein, and bends each node into the most favourable position based on score, or chooses a random direction otherwise.

**Regression**: Makes a number of random folds, only keeping the new shape if it has a better score than the previous shape.

**Simulated Annealing**: Similar to *Regression*, but has a (each iteration decreasing) chance that a new shape will be accepted even if it is worse.

**Depth First**: Recursively loops through the protein, and tries all directions in a depth-first approach. Tries directions in an order based on selected heuristics and prunes branches if they're not promising enough.
