l = [('DepthFirst - []', -17.75),
 ("DepthFirst - ['FoldAmount']", -17.725),
 ("DepthFirst - ['MinimiseDimensions']", -18.95),
 ("DepthFirst - ['Potential']", -17.2),
 ("DepthFirst - ['FoldAmount', 'MinimiseDimensions']", -15.825),
 ("DepthFirst - ['FoldAmount', 'Potential']", -20.525),
 ("DepthFirst - ['MinimiseDimensions', 'Potential']", -19.5),
 ("DepthFirst - ['FoldAmount', 'MinimiseDimensions', 'Potential']", -15.55),
 ('SimulatedAnnealingHeuristics - []', -8.675),
 ("SimulatedAnnealingHeuristics - ['FoldAmount']", -8.725),
 ("SimulatedAnnealingHeuristics - ['MinimiseDimensions']", -8.95),
 ("SimulatedAnnealingHeuristics - ['Potential']", -7.9),
 ("SimulatedAnnealingHeuristics - ['FoldAmount', 'MinimiseDimensions']", -8.35),
 ("SimulatedAnnealingHeuristics - ['FoldAmount', 'Potential']", -9.875),
 ("SimulatedAnnealingHeuristics - ['MinimiseDimensions', 'Potential']", -7.075),
 ("SimulatedAnnealingHeuristics - ['FoldAmount', 'MinimiseDimensions', "
  "'Potential']",
  -8.05)]

import matplotlib.pyplot as plt

# Data
data = [
    ('DepthFirst', [], -17.8),
    ('DepthFirst', ['F'], -17.7),
    ('DepthFirst', ['D'], -19.0),
    ('DepthFirst', ['P'], -17.2),
    ('DepthFirst', ['F', 'D'], -15.8),
    ('DepthFirst', ['F', 'P'], -20.5),
    ('DepthFirst', ['D', 'P'], -19.5),
    ('DepthFirst', ['F', 'D', 'P'], -15.6),
    ('SimulatedAnnealingHeuristics', [], -8.7),
    ('SimulatedAnnealingHeuristics', ['F'], -8.7),
    ('SimulatedAnnealingHeuristics', ['D'], -9.0),
    ('SimulatedAnnealingHeuristics', ['P'], -7.9),
    ('SimulatedAnnealingHeuristics', ['F', 'D'], -8.4),
    ('SimulatedAnnealingHeuristics', ['F', 'P'], -9.9),
    ('SimulatedAnnealingHeuristics', ['D', 'P'], -7.1),
    ('SimulatedAnnealingHeuristics', ['F', 'D', 'P'], -8.1)
]

depth_first_data = [item for item in data if item[0] == 'DepthFirst']
sim_anneal_data = [item for item in data if item[0] == 'SimulatedAnnealingHeuristics']


df_none_value = next((value for name, heuristic, value in depth_first_data if not heuristic), None)
sa_none_value = next((value for name, heuristic, value in sim_anneal_data if not heuristic), None)

plt.figure(figsize=(10, 6))

# Plotting for DepthFirst
for name, heuristic, value in depth_first_data:
    label = ''.join(heuristic) if heuristic else 'None'
    plt.plot(label, value, 'ro', label='DepthFirst' if label == 'None' else '') # Red for DepthFirst
    plt.text(label, value, f"{value:.1f}", ha='right')

# Plotting for SimulatedAnnealingHeuristics
for name, heuristic, value in sim_anneal_data:
    label = ''.join(heuristic) if heuristic else 'None'
    plt.plot(label, value, 'bo', label='SimulatedAnnealing' if label == 'None' else '') # Blue for SimulatedAnnealingHeuristics
    plt.text(label, value, f"{value:.1f}", ha='left')

# Adding horizontal lines for 'None' values without labels
if df_none_value is not None:
    plt.axhline(df_none_value, color='red', linestyle='dotted')

if sa_none_value is not None:
    plt.axhline(sa_none_value, color='blue', linestyle='dotted')

# Adding details to the plot
plt.axhline(0, color='black', lw=2) # Timeline
plt.title('Comparison of DepthFirst and SimulatedAnnealing using Heuristics')
plt.xlabel('Heuristics')
plt.ylabel('Average of 250 Proteins length= 16')
plt.grid(True)

# Adding a legend
plt.legend()

# Show the plot
plt.savefig("./output/DF_SA_with_Heuristics.png")
plt.show()
