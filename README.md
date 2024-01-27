# MinProg-AH: Protein Po(w)der
Jarec Schouten, Joël van Ingen & Wolf Gautier

## Algorithms
**Brute Force**: Tries all possible combinations. Given enough time, this algorithm will find the best shape, but is impractical to run for sufficiently large proteins.

**PureRandom**: Generates a random list of directions, ensures it's valid and computes the resulting score.

**IterativeRandom**: Loops through the protein, and finds the available directions for each node separately, trying to prevent collisions.

**Spiral**: Arranges the protein in a spiral shape.

**Greedy**: Loops through the protein, and bends each node into the most favourable position based on score, or chooses a random direction otherwise.

**Regression**: Makes a number of random folds, only keeping the new shape if it has a better score than the previous shape.

**Simulated Annealing**: Similar to *Regression*, but has a (each iteration decreasing) chance that a new shape will be accepted even if it is worse.

**Depth First**: Recursively loops through the protein, and tries all directions in a depth-first approach. Tries directions in an order based on selected heuristics and prunes branches if they're not promising enough.


## Assingment results
TODO: ADD PLOTS AND MAKE PRETTY

### 2 Dimensions
### 2.A
Depth-First: Found score of -6 in under 1000 iterations using alpha=0.5 and beta=6.

### 2.B
Depth-First: Found score of -8 in 5800 iterations using alpha=0.5 and beta=7, and found score of -9 in 69200 iterations using alpha=0.4 and beta=8.

### 2.C
