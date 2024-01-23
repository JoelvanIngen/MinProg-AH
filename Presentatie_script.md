## Notities

* Onthoud dat 3D werken al "impressive" is mbt deze opdracht.
* Distributies van alle algortimen
* Heuristieken

## Oefenpresentatie

### Intro (Joel/Jarec/Wolf)

* Voorstellen

#### Uitleg van de Case (Joel/Jarec)

* Benoem het HP-Proteinfolding problem (Protein Powder)

* Protein --> Nodes (Aminozuren)
* Aminozuren: Hydrofoob (H), Polaire (P), Cysteine(C)
* Leg uit dat het een minimalisatie probleem is

#### Hard- en softconstraints (Jarec)

* Hard constraint: Proteine mag zichzelf niet kruisen
* Soft constraint: Score moet zo laag mogelijk zijn
* "Soft constraint: Zo veel mogelijk (H) - (C) connecties?"

### Eerste Stappen

#### State-space en Data-structuur (Jarec)

* Datastructuur bestaat voornamelijk uit Protein, Nodes en 3DVecs. BELANGRIJK Directie ten opzichte van vorige node.
* State-space: 3^(n-1). Waarbij n het aantal aminozuren is. Dit wil niet zeggen dat alle configuraties legaal zijn.

#### Baseline (Joel)

* Begonnen met PureRandom en IterativeRandom (Laat histogrammen zien) 
* Spiral en Greedy volgde al snel (Geef voorbeeld outputs)

### Algoritmen

#### Random Regression & Simulated Annealing (Wolf)
* Vooral uitkomst gericht, mooie grafieken.
* Uiteleggen wat voor invloed de heuristieken hebben op Simulated Annealing


#### Iterative Greedy (Depth-First) (Joel)


### Heuristieken

#### FoldAmount/Potential/MinimiseDimension (Joel)

### Afsluiter

#### Animatie nu? (Wolf)




### Toekomst

#### Bruteforce (Jarec) & AmadeusFold (Wolf)


### Bonus slides

* Goede en slechte folds
* 3D Visualisatie? (Wolf)

