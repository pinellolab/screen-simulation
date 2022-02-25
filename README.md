# screen-simulation
Simulate bulk CRISPR screen with sorting schemes and reproter/target edits

### Setup  
In the cloned folder, run  
```
pip install -e .
```

Then you can import the `screen_simulation` module in python or jupyter notebook.  
```
from screen_simulation.SimulatedScreen import *  
```

### Usage
#### Generating new simulation of a screen  
```
sim = SimulatedScreen()
```

Simulating a single replicate  
```
sim.simulate_rep()
```

Simulating multiple replicates  
```
sim.simulate_reps()
```

Above is equivalent to simulating individual experimental steps:  
```
sim.simulate_cells()    # simulate cells with phenotype and target site edits
sim.sort_cells()        # sort cells
sim.get_genomes()       # get genome (only few genomes are obtained in experimental sampling for PCR)
sim.amplify_reads()     # amplify reads
sim.get_read_counts()   # get read counts
```


#### Screen information can be retrieved  
```
$sim
SimulatedScreen object with 
	sorting_mode = topbot
	sorting_bins = ['bot', 'mid', 'top']
	sample_names = ['bot', 'mid', 'top', 'bulk']
	sorting_bins_tested = ['bot', 'top']
	n_reps = 4
	n_targets = 700
	n_guides_per_target = 5
	n_guides = 3500
	edit_rate_dist = [0.5178685  0.60745907 0.45273832 ... 0.2441598  0.23563944 0.27550624]
	guide_coverage_dist = [1829.  983. 3272. ... 1081.  582.  665.]
	effect_size_dist = [0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
 0.  0.  0.  0.  0.  0.  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1. ]
	n_total_cells = 8000000
	n_cells_per_rep = 2000000
	n_bulk_cells = 2666666
	n_genome_per_sample = 100000
	n_reads_per_sample = 2000000
	has_reporter = False
	measures = ['guide', 'target_edit']
	measures_all = ['guide', 'target_edit', 'guide_UMI', 'target_edit_UMI']
```

```
$ sim.guide_info
$ sim.cells
$ sim.samples
$ sim.screen_res # result of sim.simulate_reps() is saved here
```

## Setting custom parameters for simulation
### Effect size distribution
```
$ import numpy as np
$ effect_size_dist = np.concatenate((np.linspace(-1, -0.1, 10), np.repeat(0, 50), np.linspace(0.1, 1, 10)))
$ sim = SimulatedScreen(effect_size_distribution = effect_size_dist)
$ sim
SimulatedScreen object with 
	sorting_mode = topbot
	sorting_bins = ['bot', 'mid', 'top']
	sample_names = ['bot', 'mid', 'top', 'bulk']
	sorting_bins_tested = ['bot', 'top']
	n_reps = 4
	n_targets = 700
	n_guides_per_target = 5
	n_guides = 3500
	edit_rate_dist = [0.5178685  0.60745907 0.45273832 ... 0.2441598  0.23563944 0.27550624]
	guide_coverage_dist = [1829.  983. 3272. ... 1081.  582.  665.]
	effect_size_dist = [-1.  -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1  0.   0.   0.   0.
  0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
  0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
  0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
  0.   0.   0.   0.   0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1. ]
	n_total_cells = 8000000
	n_cells_per_rep = 2000000
	n_bulk_cells = 2666666
	n_genome_per_sample = 100000
	n_reads_per_sample = 2000000
	has_reporter = False
	measures = ['guide', 'target_edit']
	measures_all = ['guide', 'target_edit', 'guide_UMI', 'target_edit_UMI']
```
