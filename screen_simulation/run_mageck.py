from typing import Tuple
import matplotlib.pyplot as plt
import os
from functools import reduce
import pandas as pd

from SimulatedScreen import (
    SimulatedScreen,
    _aggregate_counts
    )


def write_mageck_input(
    screen: SimulatedScreen,
    outfile_name: str
) -> pd.DataFrame:
    measures = ["guide", "target_edit", "reporter_edit"]
    try:
        to_merge = [x.set_index(["target_id", "guide_num"])[
        ["{}_{}".format(s, m) for s in screen.samples for m in measures]
        ].add_prefix("rep{}_".format(i)) for i, x in enumerate(screen.screen_res)]
    except KeyError:
        print(screen)
    
    count_tbl = reduce(lambda left, right: left.join(right), to_merge)
    count_tbl = count_tbl.reset_index().fillna(0)

    count_tbl.insert(0, "sgRNA", count_tbl.reset_index().target_id.astype(str) + "_" + count_tbl.guide_num.astype(str))
    count_tbl = count_tbl.rename(columns = {"target_id":"gene"})
    count_tbl = count_tbl.drop('guide_num', axis = 1).astype(int)
    count_tbl.to_csv(outfile_name, sep = "\t", index = False)
    
    return(count_tbl)


def run_mageck(
    screen: SimulatedScreen,
    simulation_rep: int,
    measure: str = "guide", 
    rerun: bool = False,
    out_dir: str = None,
    ) -> str:
    """Simulate & save the simulated counts and run mageck on it."""
    assert measure in ["guide", "target_edit", "reporter_edit"]
    screen.mageck_input_name = "sorted_reads/mageck_input_{}_{}.tsv".format(
        screen.get_id(),
        simulation_rep)

    if not os.path.exists(screen.mageck_input_name) or rerun:
        _ = write_mageck_input(screen.mageck_input_name, screen)
        
    screen.mageck_res_name_prefix = "mageck_results/{}_{}".format(
        screen.get_id,
        simulation_rep
    )

    def get_samples_input(sort_bin, measure):
        return(",".join(
        ["rep{}_{}_{}".format(i, sort_bin, measure) for i in range(nreps)]))
    
    mageck_ctrl_samples = get_samples_input("bulk", measure)

    for sorting_bin in screen.sorting_bins:
        if sorting_bin == "mid": continue
        mageck_res_name = mageck_res_name_prefix + "_{}.gene_summary.txt".format(
            sorting_bin)
        mageck_treatment_samples = get_samples_input(sorting_bin, measure)
        if not os.path.exists(mageck_res_name) or rerun:
            os.system("mageck test -k "+mageck_input_name+
                      " -t "+mageck_treatment_samples + 
                      " -c "+mageck_ctrl_samples+
                  " -n " + screen.mageck_res_name_prefix + "_{}".format(sorting_bin) + 
                      " --paired")

    


def get_mageck_sensitivity(
    edit_count_stats_file,
    sorting_mode,
    nreps,
    total_cells,
    bulk_cells,
    n_reads_per_sample,
    n_iter = 10
    ) -> Tuple[pd.DataFrame]:

    sim_sensitivities = []
    mageck_res_bins_all = []

    for r in range(n_iter):
        screen = SimulatedScreen(
            edit_count_stats_file = edit_count_stats_file,
            sorting_mode = sorting_mode,
            nreps = nreps,
            total_cells = total_cells,
            bulk_cells = total_cells / 3,
            n_reads_per_sample = n_reads_per_sample
            ):

        screen.simulate_reps()
        run_mageck(screen)

        bin_sensitivities = []
        mageck_reses = [self.target_info]
        for sorting_bin in screen.sorting_bins:
            mageck_res_name = screen.mageck_res_name_prefix + "_{}.gene_summary.txt".format(
                sorting_bin)
            try:
                mageck_res = pd.read_csv(mageck_res_name, sep = "\t")
            except:
                raise IOError("ERROR in reading {}".format(mageck_res_name))
                return((None, None))
            mageck_res = screen.target_info.merge(
                mageck_res, 
                right_on = "id", 
                left_on = "target_id"
                )
            if sorting_bin in ["top", "high"]:
                mageck_res["significant"] = mageck_res["pos|fdr"] < 0.05
            elif sorting_bin in ["bot", "low"]:
                mageck_res["significant"] = mageck_res["neg|fdr"] < 0.05
            else:
                raise ValueError("Unknown sorting bin")
            bin_sensitivity = mageck_res.groupby("effect_size").agg({"significant":"mean"})
            mageck_res = mageck_res[["significant"]].rename(columns = {"significant":sorting_bin})
            mageck_reses.append(mageck_res)
            bin_sensitivities.append(bin_sensitivity)

        mageck_res_bins = pd.concat(mageck_reses, axis = 1)
        mageck_res_bins["any_bin"] = mageck_res_bins[self.sorting_bins].apply(any, axis = 1)
        mean_sens_any_bins = mageck_res_bins.groupby("effect_size").agg({"any_bin":"mean"})

        sim_sensitivity = pd.concat(bin_sensitivities + [mean_sens_any_bins],axis = 1).astype(float)
        sim_sensitivity.columns = self.sorting_bins + ["any"]
    
        sim_sensitivity["sim_iter"] = r
        sim_sensitivities.append(sim_sensitivity)
        mageck_res_bins["sim_iter"] = r
        mageck_res_bins_all.append(mageck_res_bins)

    sim_sens_df = pd.concat(sim_sensitivities, axis = 0)
    mageck_res_df = pd.concat(mageck_res_bins_all, axis = 0)
    mean_sens_df = sim_sens_df.groupby("effect_size").agg("mean")[sorting_bins + ["any"]]
    
    return((mean_sens_df, mageck_res_df))
   