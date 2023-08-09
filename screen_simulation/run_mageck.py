from typing import Tuple
import os
import subprocess
import pandas as pd
from tqdm import tqdm

from screen_simulation.SimulatedScreen import SimulatedScreen


def write_mageck_input(screen: SimulatedScreen, outfile_name: str) -> pd.DataFrame:
    """
    Write the saved count information of screen into the format of MAGeCK input.
    """
    # try:
    #     to_merge = [
    #         x.guides.set_index(["target_id", "guide_num"])[
    #             [f"{s}_{m}" for s in screen.samples for m in screen.measures_all]
    #         ].add_prefix(f"rep{i}_")
    #         for i, x in enumerate(screen.screen_res)
    #     ]
    # except KeyError as e:
    #     print("Can't read screen information from the input screen object:")
    #     print(screen)
    #     raise e

    # count_tbl = reduce(lambda left, right: left.join(right), to_merge)
    # count_tbl = count_tbl.reset_index().fillna(0)

    # count_tbl.insert(
    #     0,
    #     "sgRNA",
    #     count_tbl.reset_index().target_id.astype(str)
    #     + "_"
    #     + count_tbl.guide_num.astype(str),
    # )
    # count_tbl = count_tbl.rename(columns={"target_id": "gene"})
    # count_tbl = count_tbl.drop("guide_num", axis=1).astype(int)
    count_tbl = screen.screen_res[0].to_mageck_input()
    count_tbl.to_csv(outfile_name, sep="\t", index=False)

    print(f"MaGECK input written at {outfile_name}")
    return count_tbl


def run_mageck(
    screen: SimulatedScreen,
    simulation_rep: int,
    measure: str = "guide",
    rerun: bool = False,
    out_dir: str = None,
) -> str:
    """
    Simulate & save the simulated counts and run mageck on it.
    The output is saved under mageck_results/
    """
    assert measure == "guide", measure
    screen.mageck_input_name = "sorted_reads/mageck_input_{}_{}.tsv".format(
        screen.get_id(), simulation_rep
    )

    if rerun or not os.path.exists(screen.mageck_input_name):
        _ = write_mageck_input(screen, screen.mageck_input_name)

    screen.mageck_res_name_prefix = "mageck_results/{}_{}".format(
        screen.get_id(), simulation_rep
    )

    def get_samples_input(sort_bin, measure):
        return ",".join([f"rep{i}_{sort_bin}" for i in range(screen.n_reps)])

    mageck_ctrl_samples = get_samples_input("bulk", measure)

    for sorting_bin in screen.sorting_bins:
        mageck_res_name = (
            screen.mageck_res_name_prefix
            + "_{}_{}.gene_summary.txt".format(sorting_bin, measure)
        )
        mageck_treatment_samples = get_samples_input(sorting_bin, measure)
        if not os.path.exists(mageck_res_name) or rerun:
            result = subprocess.run(
                f"source activate jy_anbe_py38; mageck test -k {screen.mageck_input_name} -t {mageck_treatment_samples} -c {mageck_ctrl_samples} -n {screen.mageck_res_name_prefix}_{sorting_bin}_{measure} --paired",
                capture_output=True,
                shell=True,
            )
            if result.returncode != 0:
                print(result.stderr)


def get_mageck_sensitivity(
    *args, FDR_cutoff=0.05, n_iter=10, rerun=False, **kwargs
) -> Tuple[pd.DataFrame]:
    """
    Get average sensitivity for each effect size from n_iter simulations.
    For a given condiiton, this method will calculate
    """
    sim_sensitivities = []
    mageck_res_bins_all = []

    for r in tqdm(range(n_iter)):
        screen = SimulatedScreen(*args, **kwargs)

        screen.simulate_reps()

        for measure in ["guide"]:
            bin_sensitivities = []
            mageck_reses = [screen.target_info]
            run_mageck(screen, r, measure=measure, rerun=rerun)
            for sorting_bin in screen.sorting_bins_tested:
                mageck_res_name = (
                    screen.mageck_res_name_prefix
                    + f"_{sorting_bin}_{measure}.gene_summary.txt"
                )
                try:
                    mageck_res = pd.read_csv(mageck_res_name, sep="\t")
                except IOError:
                    raise IOError(f"ERROR in reading {mageck_res_name}")
                mageck_res = screen.target_info.merge(
                    mageck_res, right_on="id", left_on="target_id"
                ).fillna(False)
                if sorting_bin in ["top", "high"]:
                    mageck_res["significant"] = mageck_res["pos|fdr"] < FDR_cutoff
                elif sorting_bin in ["bot", "low"]:
                    mageck_res["significant"] = mageck_res["neg|fdr"] < FDR_cutoff
                else:
                    raise ValueError("Unknown sorting bin")
                bin_sensitivity = mageck_res.groupby("effect_size").agg(
                    {"significant": "mean"}
                )
                mageck_res = mageck_res.set_index("id", drop=True)[
                    ["significant"]
                ].rename(columns={"significant": sorting_bin})
                mageck_reses.append(mageck_res)
                bin_sensitivities.append(bin_sensitivity)

            mageck_res_bins = pd.concat(mageck_reses, axis=1)  # 700 x 5
            mageck_res_bins["any_bin"] = mageck_res_bins[
                screen.sorting_bins_tested
            ].apply(
                any, axis=1
            )  # 700 x 1
            mean_sens_any_bins = mageck_res_bins.groupby("effect_size").agg(
                {"any_bin": "mean"}
            )  # 11
            assert len(mean_sens_any_bins) == len(
                bin_sensitivities[0]
            ), bin_sensitivities[0]
            sim_sensitivity = pd.concat(
                bin_sensitivities + [mean_sens_any_bins], axis=1
            ).astype(float)
            try:
                sim_sensitivity.columns = screen.sorting_bins_tested + ["any"]
            except:
                print(sim_sensitivity)
                print(screen)

            sim_sensitivity["sim_iter"] = r
            sim_sensitivity["measure"] = measure
            sim_sensitivities.append(sim_sensitivity)
            mageck_res_bins["sim_iter"] = r
            mageck_res_bins["measure"] = measure
            mageck_res_bins_all.append(mageck_res_bins)

    sim_sens_df = pd.concat(sim_sensitivities, axis=0)
    mageck_res_df = pd.concat(mageck_res_bins_all, axis=0)
    mean_sens_df = sim_sens_df.groupby(["measure", "effect_size"]).agg("mean")[
        screen.sorting_bins_tested + ["any"]
    ]

    return (mean_sens_df, mageck_res_df)
