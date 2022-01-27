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
sim
#SimulatedScreen object with 
#	sorting_mode = topbot
#	sorting_bins = ['bot', 'mid', 'top']
#	sample_names = ['bot', 'mid', 'top', 'bulk']
#	n_reps = 4
#	n_targets = 700
#	n_guides_per_target = 5
#	n_guides = 3500
#	n_total_cells = 8000000
#	n_cells_per_rep = 2000000
#	n_bulk_cells = 2666666
#	n_genome_per_sample = 100000
#	n_reads_per_sample = 2000000
#	vars_per_mu = 10
#	mu_steps = 10
#	max_effect_size = 1
#	has_reporter = False
```

```
sim.guide_info
sim.cells
sim.samples
sim.screen_res # result of sim.simulate_reps() is saved here
```
