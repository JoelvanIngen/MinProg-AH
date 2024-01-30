# Assignments
**PureRandom**: Best of N runs
**IterativeRandom**: Best of N runs
**Spiral**: Single run (deterministic)
**Greedy**: Best of 1000 runs (semi-deterministic)
**SimulatedAnnealing**: Best after running for N seconds

## 2 - 2D
### A: HHPHHHPHPHHHPH
**PureRandom**: -6 (Best of 100000 runs)
**IterativeRandom**: -6 (Best of 50000 runs)
**Spiral**: -4
**Greedy**: -6
**DepthFirst**: -6, Alpha=0.5, Beta=6
**SimulatedAnnealing** -6

### B: HPHPPHHPHPPHPHHPPHPH
**PureRandom**: -7 (Best of 100000 runs)
**IterativeRandom**: -8 (Best of 50000 runs)
**Spiral**: -4
**Greedy**: -7
**DepthFirst**: -8, Alpha=0.5, Beta=7
**SimulatedAnnealing**: -8

### C: PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP
**PureRandom**: -7 (Best of 100000 runs)
**IterativeRandom**: -10 (Best of 50000 runs)
**Spiral**: -5
**Greedy**: -11
**DepthFirst**: -10, Alpha=0.5, Beta=7
**SimulatedAnnealing**: -9

### D: HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH
**PureRandom**: -10 (Best of 10000 runs)
**IterativeRandom**: -11 (Best of 10000 runs)
**Spiral**: -10
**Greedy**: -15
**DepthFirst**: -16, Alpha=0.5, Beta=11
**SimulatedAnnealing**: -22 

## 3 - 3D
### A: PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP
**PureRandom**: -14 (Best of 50000 runs)
**IterativeRandom**: -18 (Best of 50000 runs)
**Spiral**: -10
**Greedy**: -19
**DepthFirst**: -19, Alpha=0.45, Beta=12 (100000 iterations)
**SimulatedAnnealing**: -25

### B: CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC
**PureRandom**: -27 (Best of 50000 runs)
**IterativeRandom**: -29 (Best of 50000 runs)
**Spiral**: -16
**Greedy**: -31
**DepthFirst**: -19, Alpha=0.6, Beta=20 (100000 iterations)
**SimulatedAnnealing**: -44

### C: HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH
**PureRandom**: -16 (Best of 10000 runs)
**IterativeRandom**: -23 (Best of 10000 runs)
**Spiral**: -12
**Greedy**: -25
**DepthFirst**: -25, Alpha=0.6, Beta=25 (100000 iterations)
**SimulatedAnnealing**: -36

### D: HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH
**PureRandom**: -15 (Best of 10000 runs)
**IterativeRandom**: -21 (Best of 10000 runs)
**Spiral**: -15
**Greedy**: -25
**DepthFirst**: -29, Alpha=0.5, Beta=15 (100000 iterations)
**SimulatedAnnealing**: -31

## 4 - Analyse
Zoals we in het begin van het project voorspelden schaalt het verschil tussen 
de baselines en de algoritmen met de state space van de sequence. Dit 
blijkt ten eerste uit dat het verschil tussen de scores van de baseline en van 
sommige algoritmen toeneemt naarmate de sequence length toeneemt, maar ook uit 
het nog grotere verschil tussen de baselines en de algoritmen in 3 dimensies. 
Verder komen ook voor- en nadelen van algoritmes naar voren, die belangrijk 
worden als je folds zou willen renderen met een specifiek doel; het depth-first 
algoritme is rigoureuzer en levert waarschijnlijk als je onbeperkte runtime 
hebt een lagere score op, maar is dusdanig veel langzamer dan simulated 
annealing dat het op sommige plekken in de bovenstaande resultaten ingehaald 
wordt binnen de beperkte rekentijd die beide algoritmes kregen. Voor een 
toepassing waar de score het belangrijkst is zou dus een rigoureuzer algoritme 
als depth-first een betere keuze zijn, terwijl andere opties als simulated 
annealing juist voordeliger zouden kunnen zijn als quantiteit en snelheid 
belangrijker is.
