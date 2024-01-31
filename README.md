
# MinProg-AH: Protein Po(w)der
Jarec Schouten, Joël van Ingen & Wolf Gautier

<div align="center">
<figure>
    <img src="./trophies/atest-cropped2.gif" height="439">
    <h4></h4>
    <figcaption>A sample 3D protein fold.</figcaption>
</figure>
</div>

## Table of Contents

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [Case Explanation](#case-explanation)
- [Running the program](#running-the-program)
- [Configurations](#configurations)
   * [Algorithms](#algorithms)
   * [Parameters](#parameters)
   * [Heuristics: Only available in the DepthFirst algorithm](#heuristics-only-available-in-the-depthfirst-algorithm)
- [ML approach: AmadeusFold](#amadeusfold)
   * [Pipeline overview](#pipeline-overview)
- [Example Run and Output](#example-run-and-output)
- [Separate Experiments (milestone)](#separate-experiments-milestone)

<!-- TOC end -->


<!-- TOC --><a name="case-explanation"></a>
## Case Explanation
In our chosen case, we implement solutions for efficiently solving the H-P protein model. 
This a simplified model for folding protein structures in order to achieve a low energy state i.e., form the most H-bonds between atoms in the protein.
These bonds can only form if the modecules are directly next to eachother on a grid.
However, getting the right molecules positioned next to eachother requires a well-structured protein, as the protein cannot cross itself.

We implemented multiple algorithms that are designed to achieve these low-energy folds in a short timespan, and we extended the model to 3D.

<!-- TOC --><a name="running-the-program"></a>
## Running the program
By running main.py, there is a small UI in the terminal that allows the user to enter a sequence and select an algorithm alongside paramaters and heuristcs.
More information on these configurations is read below.

<!-- TOC --><a name="configurations"></a>
## Configurations

<!-- TOC --><a name="algorithms"></a>
### Algorithms
- **Brute Force**: Tries all possible combinations. Given enough time, this algorithm will find the best shape, but is impractical to run for sufficiently large proteins.

- **PureRandom**: Generates a random list of directions, ensures it's valid and computes the resulting score.

- **IterativeRandom**: Loops through the protein, and finds the available directions for each node separately, trying to prevent collisions.

- **Spiral**: Arranges the protein in a spiral shape.

- **Greedy**: Loops through the protein, and bends each node into the most favourable position based on score, or chooses a random direction otherwise.

- **Regression**: Makes a number of random folds, only keeping the new shape if it has a better score than the previous shape.

- **Simulated Annealing**: Similar to *Regression*, but has a (each iteration decreasing) chance that a new shape will be accepted even if it is worse.

- **Depth First**: Recursively loops through the protein, and tries all directions in a depth-first approach. Tries directions in an order based on selected heuristics and prunes branches if they're not promising enough.

- **Beam Search**: Explores position that has the highest points at that moment. Saves others in a queue.

<!-- TOC --><a name="parameters"></a>
### Parameters

- **Dimensions**: The user can choose whether they want to run the algorithm in 2D or 3D.

- **Debugging**: Outputs internal data describing at where the algorithm is at that point. (ONLY USE IF ALGORITHM CRASHES)

- **Verbose**: Outputs more global information, i.e. when the algorithm has found a new best score.

- **Progressbar**: Displays a progressbar while running.

- **State History**: Allows the algorithm to save previous states to enable the user to animate the algorithm.

- **Max Iterations**: Caps the amount of iterations the algorithm is allowed to run for.

<!-- TOC --><a name="heuristics-only-available-in-the-depthfirst-algorithm"></a>
### Heuristics: Only available in the DepthFirst algorithm

- **MinimiseDimensions**: Nudges the algorithm in directions where the overall dimensions of the protein are kept as small as possible.

- **PotentialPlus**: Nudges the algorithm in directions where hydrophobic nodes are kept together and polar nodes are placed further away.

- **FoldAmount**: Nudges the algorithm in directions where the protein makes the most amount of 90 degree folds.

<!-- TOC --><a name="amadeusfold"></a>
## ML approach: AmadeusFold
In addition to the mentioned algorithms, an LSTM-based machine learning approach was implemented. Once trained, it takes a sequence as input and outputs an order, following the conventions of the rest of the project. 

<!-- TOC --><a name="pipeline"></a>
### Pipeline overview
In general, the following steps should be taken to obtain a functioning model.
1. **Generate a dataset:** using the file data_generator.py, a customisable dataset can be generated consisting of sequences and their optimal solutions, obtained through the bruteforce algorithm. It is stored as a .csv.
2. **Choose a model structure:** models can be trained with standard loss functions, resulting in quicker training but slightly lower accuracy, or a custom loss function, yielding higher accuracy but taking longer to train. A further choice can be made between opting for a classifier architecture or one that outputs coordinates directly.
3. **Train the chosen model with the generated dataset:** once data is generated and a model is chosen, adjusting the desired parameters in a training_(model).py script begins training the model and saves the state dict once training finishes.
4. **Obtaining fold results:** within either run.py or run_classifier.py, simply adjust the sequence variable to the one that should be folded and run the script. An order or a set of coordinates will be displayed to the screen.


<!-- TOC --><a name="example-run-and-output"></a>
## Example Run and Output

Running `main.py` looks as followed:
```
$ python3 main.py 
Enter sequence: PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP
Available algorithms:
[0] PureRandom
[1] IterativeRandom
[2] Spiral
[3] BruteForce
[4] Regression
[5] SimulatedAnnealing
[6] Greedy
[7] BeamSearch
[8] DepthFirst
Choose an algorithm 0-8: 8
Choose a number of dimensions (2 / 3): 2
Change parameter max_iterations (default 5000)? (y/n): y
Set parameter max_iterations to: 20000
Change parameter prune_alpha (default 0.5)? (y/n): n
Change parameter prune_beta (default 15)? (y/n): y
Set parameter prune_beta to: 12
Use debugging? (y/n): n
Keep state history for animating? (y/n): n
Use progressbar? (y/n): y
Output verbose information? (y/n): n
Add heuristic MinimiseDimensions? (y/n): y
Add heuristic PotentialPlus? (y/n): n
Add heuristic FoldAmount? (y/n): n
 98%|█████████▊| 19691/20000 [00:04<00:00, 4562.28it/s]

Final score: -19
20001it [00:05, 3966.01it/s]                           

Process finished with exit code 0
```

Once this is run, the user will be able to find a plot of the protein in the `./experiments/output` folder.
This particular sequence and configuration yields the following protein:

<div align="center">
<figure>
    <img src="./trophies/presentationPlots/evaluate_DepthFirst_len36_dim2_a0.5_b12.png" width="640" height="480">
    <h4></h4>
    <figcaption>Output DepthFirst Algorithm</figcaption>
</figure>
</div>

<!-- TOC --><a name="separate-experiments-milestone"></a>
## Separate Experiments (milestone)
To run timed tests on specific algorithms (i.e. simulated annealing), run experiments/iterate_sa.py, adjust sequence, timing and parameters as necessary.
