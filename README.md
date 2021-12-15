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

#### Screen information can be retrieved
```
sim.guide_info
sim.cells
sim.samples
sim.screen_res # result of sim.simulate_reps() is saved here
```
