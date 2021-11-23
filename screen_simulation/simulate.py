from multiprocessing import Pool 
from SimulatedScreen import SimulatedScreen
from run_mageck import run_mageck



def simulate_screens(
    
    ):
    
    reses = []
    with Pool(100) as pool:
    iterable = [(i, 'x') for i in range(10)]
    for _ in tqdm(pool.istarmap(get_mageck_sensitivity, zipped_params),
                       total=len(zipped_params)):
        reses.append(_)

