from typing import List, Optional, Union
from typing import Optional, Union, List
from inspect import ismethod

import numpy as np
import torch
import pandas as pd
from pprint import pprint
import scipy.stats
import bean as be
from bean import ReporterScreen
from bean.model.utils import scale_pi_by_accessibility


def _aggregate_counts(df, prefix=None, group_by=["target_id", "guide_num"]):
    counts = df.groupby(group_by).agg(
        {"edited": ["count", "sum"], "reporter_edited": "sum"}
    )

    counts.columns = ["guide", "target_edit", "reporter_edit"]
    if not prefix is None:
        counts.columns = counts.columns.map(lambda s: "{}_".format(prefix) + s)
    return counts


class SimulatedScreen:
    """
    Simulated screen with phenotype variances and sorting.
    
    Attributes:
    edit_count_stats_file (str): .csv file with name, edit_rate, and guide_count columns. The edit rate and guide counts will be sampled from this dataframe. \
        Assumed to have the same number of total guides as the n_guides = n_targets*n_guides_per_target.
    sorting_mode (str): Sorting scheme: "topbot" - top 30 / bot 30 , "bins" - 20,40,60,80,100 quantiles
    nreps (int): Number of replicates that will be used (each replicate assumed to be iid)
    n_targets (int): Number of target genes/variants.
    n_guides_per_target (int): Number of guides per target.
    n_total_cells (int): Total number of cells that will be sorted.
    n_bulk_cells (int): Number of cells in the bulk (control, unsorted) sample.
    n_reads_per_sample (int): Guide read depth of each sample.
    vars_per_mu (int) Number of target gene/variants with the same effect size (mu).
    mu_steps (int) Number of steps that will be used
    max_effect_size (float) -- Largest effect size (z-score) of phenotype of target variants in the screen.
    screen_res (List[pd.DataFrame])-- List of read count tables of sorted/bulk samples.
    scale_by_accessibility -- Scale reporter editing rate by accessibility to obatin endogenous target site. When np.array is provided, use the value as the raw accessibility score. When True is passed, accessibility from data is provided.
    scale_with_variability -- scale endogeneous site edit rate by reporter edit rate with variability observed from data
    """

    edit_rates_by_target_saved_data = np.loadtxt(
        "/data/pinello/PROJECTS/2021_08_ANBE/data/all_ldlvar_edit_rates_by_target.csv"
    )
    edit_rates_saved_data = np.loadtxt(
        "/data/pinello/PROJECTS/2021_08_ANBE/data/all_ldlvar_edit_rates.csv"
    )
    guide_coverage_saved_data = np.loadtxt(
        "/data/pinello/PROJECTS/2021_08_ANBE/data/all_ldlvar_bulk_coverage_sums.csv"
    )
    accessibility_saved_data = np.loadtxt(
        "/data/pinello/PROJECTS/2021_08_ANBE/data/gRNA_info/20221215_LDLvar_target_accessibility.csv",
        delimiter=",",
    )  # (n_targets,)
    endo_by_reporter_residual_saved_data = np.loadtxt(
        "/data/pinello/PROJECTS/2021_08_ANBE/data/accessibility/20221215_ratio_residual.csv"
    )  # scale endo/reporter by this factor

    def __init__(
        self,
        n_targets: int = 700,
        n_guides_per_target: int = 5,
        n_total_cells: int = 10**6 * 8,
        n_bulk_cells=None,
        selection: str = "sorting",
        sorting_mode: str = "topbot",
        survival_sample_times: Optional[List[float]] = None,
        nreps: int = 4,
        n_genome_per_sample: int = int(10e4),
        n_reads_per_sample: int = 10**6 * 2,
        effect_size_distribution: Optional[Union[np.ndarray, str]] = None,
        guide_coverage_distribution: Optional[Union[np.ndarray, str]] = "data",
        edit_rate_distribution: Optional[Union[np.ndarray, str]] = "data",
        has_reporter=False,
        sample_edit_rate_by_target: bool = True,
        scale_by_accessibility: Optional[Union[np.ndarray, bool]] = False,
        scale_with_variability: bool = False,
    ):
        if selection == "sorting":
            self.sorting_mode = sorting_mode
            if self.sorting_mode == "topbot":
                self.sorting_bins = ["bot", "mid", "top"]
            elif self.sorting_mode == "bins":
                self.sorting_bins = ["bot", "low", "mid", "high", "top"]
            else:
                raise ValueError("Invalid sorting_mode: {}".format(sorting_mode))
        elif selection == "survival":
            self.survival_sample_times = np.array(survival_sample_times)
        else:
            raise ValueError("Invalid selection: {}".format(selection))
        self.sample_names = self.sorting_bins + ["bulk"]
        self.sorting_bins_tested = [sb for sb in self.sorting_bins if sb != "mid"]

        self.n_reps = nreps
        self.n_targets = n_targets
        self.n_guides_per_target = n_guides_per_target
        self.n_guides = self.n_targets * self.n_guides_per_target
        self.sample_edit_rate_by_target = sample_edit_rate_by_target
        self.scale_with_variability = scale_with_variability
        self.scale_by_accessibility = scale_by_accessibility
        if edit_rate_distribution is None:
            self.edit_rate_distribution = np.array([1])

        elif edit_rate_distribution == "data":
            if sample_edit_rate_by_target:
                self.edit_rate_distribution = type(self).edit_rates_by_target_saved_data
            else:
                self.edit_rate_distribution = type(self).edit_rates_saved_data
        elif isinstance(edit_rate_distribution, np.ndarray):
            self.edit_rate_distribution = edit_rate_distribution
        else:
            raise ValueError(
                "Invalid edit_rate_distribution: {}".format(edit_rate_distribution)
            )

        if guide_coverage_distribution is None:
            self.guide_coverage_dist = np.array([1])
        elif guide_coverage_distribution == "data":
            self.guide_coverage_dist = type(self).guide_coverage_saved_data
        elif isinstance(guide_coverage_distribution, np.ndarray):
            self.guide_coverage_dist = guide_coverage_distribution
        else:
            raise ValueError(
                "Invalid guide_coverage_distribution: {}".format(
                    guide_coverage_distribution
                )
            )

        if effect_size_distribution is None:
            self.effect_size_distribution = np.concatenate(
                (np.repeat(0, 60), np.linspace(0.1, 1, 10))
            )
        elif isinstance(effect_size_distribution, np.ndarray):
            self.effect_size_distribution = effect_size_distribution
        else:
            raise ValueError(
                "Invalid effect_size_distribution: {}".format(effect_size_distribution)
            )

        # read accessibility and get scale
        if isinstance(scale_by_accessibility, np.ndarray):
            self.target_accessibility = scale_by_accessibility
        elif scale_by_accessibility:
            self.target_accessibility = type(self).accessibility_saved_data
            self.edit_rate_distribution = type(self).edit_rates_by_target_saved_data[
                20:
            ]  # remove control
        else:
            self.target_accessibility = None

        self.n_total_cells = n_total_cells
        self.n_cells_per_rep = int(self.n_total_cells / self.n_reps)
        if n_bulk_cells is None:
            self.n_bulk_cells = int(self.n_total_cells / 3)
        else:
            self.n_bulk_cells = n_bulk_cells
        self.n_genome_per_sample = n_genome_per_sample
        self.n_reads_per_sample = n_reads_per_sample

        self.has_reporter = has_reporter
        self.measures = ["guide", "target_edit"]
        if self.has_reporter:
            self.measures.append("reporter_edit")
        measures_umi = list(map(lambda s: s + "_UMI", self.measures))
        self.measures_all = self.measures + measures_umi
        self.screen_res = []

        self.covs = None
        self.edit_rates = None
        self.effect_sizes = None

    def __repr__(self):
        out_str = "SimulatedScreen object with \n"
        for k, v in getattr(self, "__dict__").items():
            if k.startswith("_"):
                continue
            if ismethod(v):
                continue
            if isinstance(v, pd.DataFrame):
                continue
            if k in ["samples", "screen_res"]:
                continue
            out_str += "\t{} = {}\n".format(k, v)
        return out_str

    def get_id(self):
        params = "{}_{}_{}_{}".format(
            self.sorting_mode,
            self.n_cells_per_rep,
            self.n_reps,
            self.n_reads_per_sample,
        )
        return params

    def _sample_all_from_data(self, sample_by_target=True):
        """Samples (guide_count, (reporter) edit_rate, endo_edit_rate) from self.count_edit_stats"""
        if sample_by_target:
            if self.scale_by_accessibility:
                er_target_idx = np.random.choice(
                    len(self.target_accessibility), size=self.n_targets, replace=True
                )
            else:
                er_target_idx = np.random.choice(
                    self.edit_rate_distribution.shape[0],
                    size=self.n_targets,
                    replace=True,
                )
            edit_rate = self.edit_rate_distribution[er_target_idx, :].flatten()
        else:
            edit_rate = np.random.choice(
                self.edit_rate_distribution, size=self.n_guides, replace=True
            )
        guide_count_raw = np.random.choice(
            self.guide_coverage_dist, size=self.n_guides, replace=True
        )
        guide_count_norm = np.floor(
            guide_count_raw / self.guide_coverage_dist.sum() * self.n_total_cells
        )

        if not self.target_accessibility is None:
            target_acc = self.target_accessibility[er_target_idx]
            target_acc[np.isnan(target_acc)] = np.nanmedian(
                self.target_accessibility
            )  # TODO revise later?
            acc = np.repeat(target_acc, self.n_guides_per_target)
            assert acc.shape == (self.n_guides,), (
                acc.shape,
                er_target_idx.shape,
                self.target_accessibility[er_target_idx].shape,
            )
            er = torch.as_tensor(edit_rate)
            allele_rate = (
                torch.stack([1 - er, er], axis=-1)
                .unsqueeze(0)
                .unsqueeze(0)
                .expand(self.n_reps, 1, -1, -1)
            )
            endo_edit_rate = (
                scale_pi_by_accessibility(allele_rate, torch.as_tensor(acc))
                .detach()
                .numpy()
            )[0, 0, :, -1]
        else:
            endo_edit_rate = edit_rate
            acc = None
        return (guide_count_norm, edit_rate, endo_edit_rate, acc)

    def _get_effect_sizes(self):
        effect_sizes = np.resize(self.effect_size_distribution, self.n_targets)

        self.guide_info = pd.DataFrame(
            {
                "target_id": np.repeat(
                    list(range(self.n_targets)), self.n_guides_per_target
                ),  # 0,0,0,0,0,1,1,1,1,1,2,...
                "guide_num": np.tile(
                    list(range(self.n_guides_per_target)), self.n_targets
                ),
                "effect_size": np.repeat(effect_sizes, self.n_guides_per_target),
            }
        )
        self.guide_info["guide_id"] = (
            self.guide_info.target_id.astype(str)
            + "_"
            + self.guide_info.guide_num.astype(str)
        )
        self.guide_info = self.guide_info.set_index("guide_id", drop=True)

        self.target_info = (
            self.guide_info[["target_id", "effect_size"]]
            .drop_duplicates()
            .reset_index()
        )
        return effect_sizes

    def simulate_cells(self):
        """Simulates cells as dataframe as self.cells."""
        if self.effect_sizes is None:
            # Simulate the effect sizes of the variants.
            self.effect_sizes = self._get_effect_sizes()
        if self.covs is None and self.edit_rates is None:
            # Get coverage and edit rate pairs for each guide
            (
                self.covs,
                self.edit_rates,
                self.endo_edit_rates,
                self.acc,
            ) = self._sample_all_from_data(self.sample_edit_rate_by_target)

            # Cell coverage and edit rate per guide
            self.guide_info["coverage"] = self.covs.astype(int)
            self.guide_info["edit_rate"] = self.edit_rates
            self.guide_info["endo_edit_rate"] = self.endo_edit_rates
            self.guide_info["acc"] = self.acc

        cells = self.guide_info.reindex(
            self.guide_info.index.repeat(self.guide_info.coverage)
        ).reset_index()
        cells["edited"] = scipy.stats.bernoulli.rvs(cells.endo_edit_rate)
        cells["phenotype_mean"] = cells["edited"] * cells["effect_size"]
        cells["phenotype"] = scipy.stats.norm.rvs(loc=cells.phenotype_mean, scale=1)
        # TODO: for scRNA readout this line should be different @Lucas

        if self.has_reporter:
            cells["reporter_edited"] = scipy.stats.bernoulli.rvs(cells.edit_rate)
        self.cells = cells

    def expand_cells(self):
        self.cells = NotImplemented  # TODO

    def select_cells(self):
        if self.selection == "sorting":
            self.sort_cells()
        else:
            self.prolif_cells()

    def sort_cells(self):
        """Sort cells into designated sorting scheme: bulk, 1/2/3/4"""
        if self.sorting_mode == "bins":
            self.sort_cells_quantiles(
                q=[0, 0.2, 0.4, 0.6, 0.8, 1], bin_labels=self.sorting_bins
            )
        elif self.sorting_mode == "topbot":
            self.sort_cells_quantiles(q=[0, 0.3, 0.7, 1], bin_labels=self.sorting_bins)

    def sort_cells_quantiles(self, q=[0, 0.2, 0.4, 0.6, 0.8, 1], bin_labels=None):
        """Sort cells (phenotypes) into quantile bins.

        Parameters:
        q (int or list-like of float): q Parameter for pandas.qcut. \
            Number of quantiles. 10 for deciles, 4 for quartiles, etc. \
            Alternately array of quantiles, e.g. [0, .25, .5, .75, 1.] for quartiles.
        """
        if bin_labels is None:
            bin_labels = np.arange(len(q) - 1).astype(str)
        else:
            if len(bin_labels) != len(q) - 1:
                raise ValueError(
                    "provided length of bin_labels \
                    doesn't match the number of bins made by quantile"
                )

        self.cells["sorted_bin"] = pd.qcut(self.cells.phenotype, q=q, labels=bin_labels)

    def prolif_cells(self):
        """
        t: Time of selection
        """
        for t in self.survival_sample_times:
            self.cells[f"prolif_{t}"] = self.cells.phenotype.map(
                lambda x: np.exp(x * t)
            )

    def get_genomes(self):
        """
        Samples the genomes for each sequencing samples prepped.
        Number of genomes sampled is sampled as the multinomial
        distribution with p proportional to cell counts per guide.

        Assigns self.samples attribute which is Dict[str -> pd.DataFrame]
        that stores sampled genomes per sample with UMI
        """
        self.samples = {}  # Dict[str -> pd.DataFrame]
        cell_record_to_keep = ["target_id", "guide_num", "edited"]
        if self.has_reporter:
            cell_record_to_keep.append("reporter_edited")

        # sorted bins
        if self.selection == "sorting":
            for sorting_bin, df in self.cells.groupby("sorted_bin", as_index=False):
                if self.n_genome_per_sample > len(df):
                    # All genomes are sampled
                    self.samples[sorting_bin] = df[cell_record_to_keep]
                else:
                    self.samples[sorting_bin] = df[cell_record_to_keep].sample(
                        n=self.n_genome_per_sample, replace=False
                    )

            # bulk sample
            self.samples["bulk"] = self.cells[cell_record_to_keep].sample(
                n=self.n_genome_per_sample, replace=False
            )
        # proliferation screen
        else:
            for t in self.survival_sample_times:
                self.samples[f"prolif_{t}"] = pd.DataFrame(
                    np.repeat(
                        self.cells[cell_record_to_keep],
                        self.cells[f"prolif_{t}"].values,
                        axis=0,
                    ),
                    columns=cell_record_to_keep,
                )
        # Assign UMI
        for sample_name in self.sample_names:
            self.samples[sample_name]["UMI"] = list(
                range(len(self.samples[sample_name]))
            )

    def amplify_reads(self, dispersion_fn=None):
        """
        dispersion_fn: pair of floats (c0, c1), disp = c0 + c1/q
        """
        for sample_name, df in self.samples.items():
            amplified_counts = np.ones(
                len(df)
            )  # TODO: different amplification per molecule (UMI)
            amplified_counts *= self.n_reads_per_sample / amplified_counts.sum()
            # introduce NB noise
            if dispersion_fn is None:
                dispersions = 0.01
            else:
                dispersions = dispersion_fn[0] + dispersion_fn[1] / (
                    amplified_counts + 1
                )
            amplified_counts = np.random.negative_binomial(
                n=1 / dispersions,
                p=amplified_counts / (1.0 / dispersions + amplified_counts),
            ).astype(float)

            amplified_counts *= self.n_reads_per_sample / amplified_counts.sum()
            self.samples[sample_name] = pd.DataFrame(
                np.repeat(df.values, amplified_counts.astype(int), axis=0)
            )
            self.samples[sample_name].columns = df.columns

    def get_read_counts(self, rep=0):
        agg_fn = {"edited": ["count", "sum"]}
        if self.has_reporter:
            agg_fn["reporter_edited"] = "sum"
        measures_umi = list(map(lambda s: s + "_UMI", self.measures))
        cdata_list = []
        cdata_umi_list = []

        for sample_name, df in self.samples.items():
            condit_df = pd.DataFrame(columns=["replicate", "sort"])
            condit_df.loc["rep{}_{}".format(rep, sample_name)] = [rep, sample_name]
            cdata = _counts_df_to_repscreen(
                df, agg_fn, self.guide_info, condit_df, has_reporter=self.has_reporter
            )
            cdata_umi = _counts_df_to_repscreen(
                df, agg_fn, self.guide_info, condit_df, has_reporter=self.has_reporter
            )
            cdata_list.append(cdata)
            cdata_umi_list.append(cdata_umi)
        cdata_all = be.concat(cdata_list, axis=1)
        cdata_umi_all = be.concat(cdata_umi_list, axis=1)
        return (cdata_all, cdata_umi_all)

    def simulate_rep(self, rep=0):
        self.simulate_cells()
        # TODO: growing cells and adapting for viability screen
        # self.expand_cells()
        # self.select_cells()
        self.sort_cells()
        self.get_genomes()
        self.amplify_reads()
        read_counts, umi_counts = self.get_read_counts(rep)
        return (read_counts, umi_counts)

    def simulate_reps(self):
        list_cdata = []
        list_cdata_umi = []
        for rep in range(self.n_reps):
            rep_cdata, rep_cdata_umi = self.simulate_rep(rep)
            list_cdata.append(rep_cdata)
            list_cdata_umi.append(rep_cdata_umi)
        reps_cdata = be.concat(list_cdata, axis=1)
        reps_cdata_umi = be.concat(list_cdata_umi, axis=1)
        self.screen_res = (reps_cdata, reps_cdata_umi)


def _counts_df_to_repscreen(df, agg_fn, guides, condit, has_reporter=True):
    counts = df.groupby(["target_id", "guide_num"]).agg(agg_fn)
    counts = counts.reindex(
        pd.MultiIndex.from_frame(guides[["target_id", "guide_num"]]), fill_value=0
    )
    if has_reporter:
        edits = counts[[("reporter_edited", "sum")]].to_numpy()
    else:
        edits = None
    cdata = ReporterScreen(
        X=counts[[("edited", "count")]].to_numpy(),
        X_edit=edits,
        X_bcmatch=counts[[("edited", "count")]].to_numpy(),
        guides=guides,
        samples=condit,
    )
    cdata.layers["target_edited"] = counts[[("edited", "sum")]].to_numpy()
    return cdata
