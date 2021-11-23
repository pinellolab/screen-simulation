from typing import List
import numpy as np
import pandas as pd
from pprint import pprint
import scipy.stats


def _aggregate_counts(df, prefix = None, group_by = ["target_id", "guide_num"]):
    counts = df.groupby(group_by).agg(
        {"edited":["count", "sum"], "reporter_edited":"sum"})
    
    counts.columns = ["guide", "target_edit", "reporter_edit"]
    if not prefix is None:
        counts.columns = counts.columns.map(lambda s: "{}_".format(prefix) + s)
    return(counts)



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
    """

    def __init__(self, 
        edit_count_stats_file = None,
        sorting_mode: str = "topbot",
        nreps = 4,
        n_targets = 700,
        n_guides_per_target = 5,
        n_total_cells = 10**6*8,
        n_bulk_cells = None,
        n_reads_per_sample = 10**6*2,
        vars_per_mu = 10,
        mu_steps = 10,
        max_effect_size = 1
        ):

        self.sorting_mode = sorting_mode
        if self.sorting_mode == "topbot":
            self.sorting_bins = ["bot", "mid", "top"]
        elif self.sorting_mode == "bins":
            self.sorting_bins =  ["bot", "low", "mid", "high", "top"]
        else:
            raise ValueError("Invalid sorting_mode")
        self.samples = self.sorting_bins + ["bulk"]

        self.nreps = nreps
        self.n_targets = n_targets
        self.n_guides_per_target = n_guides_per_target
        self.n_guides = self.n_targets * self.n_guides_per_target

        if edit_count_stats_file is None:
            self.count_edit_stats = pd.DataFrame({
                "guide_counts" : np.ones(len(self.n_guides)),
                "edit_rate" : 1
                }
            )
        else:
            self.count_edit_stats = pd.read_csv(edit_count_stats_file)

        self.n_total_cells = n_total_cells
        self.n_cells_per_rep = int(self.n_total_cells / self.n_reps)
        if n_bulk_cells is None:
            self.n_bulk_cells = self.n_total_cells / 3
        else:
            self.n_bulk_cells = n_bulk_cells
        self.n_reads_per_sample = n_reads_per_sample
        self.vars_per_mu = vars_per_mu
        self.mu_steps = mu_steps
        self.max_effect_size = max_effect_size
        self.phenotypes = None
        self.screen_res = []


    def __repr__(self):
        out_str = "SimulatedScreen object with \n"
        for attr in dir(self):
            out_str += "\t{} = {}\n".format(attr, getattr(self, attr))
        return(out_str)

    def get_id(self):
        params = "{}_{}_{}_{}".format(
            self.sorting_mode, 
            self.n_cells_per_rep, 
            self.nreps, 
            self.n_reads_per_sample
        )
        return(params)


    def _sample_all_from_data(self):
        row = self.count_edit_stats.sample(self.n_guides, replace = True).reset_index()
        guide_count_norm = np.floor(row.guide_count / self.count_edit_stats.guide_count.sum() * self.total_cells)
        return((guide_count_norm, row.edit_rate))

    def _get_effect_sizes(self):
        no_effect_targets = self.n_targets - self.mu_steps*self.vars_per_mu
        mu_delta = np.linspace(
            self.max_effect_size/self.mu_steps, 
            self.max_effect_size, 
            num = self.mu_steps)
        coverage = np.floor(self.total_cells / self.n_guides)
        effect_sizes = np.concatenate(
            (np.zeros(no_effect_targets), 
            np.repeat(mu_delta, self.vars_per_mu)
            ))
    
        self.guide_info = pd.DataFrame({
            "target_id" : np.repeat(list(range(self.n_targets)), self.n_guides_per_target), 
            "guide_num" : np.tile(list(range(self.n_guides_per_target)), self.n_targets), 
            "effect_size" : np.repeat(effect_sizes, self.n_guides_per_target)
            })

        self.target_info = self.guide_info[["target_id", "effect_size"]].drop_duplicates()
        return(effect_sizes)


    def simulate_phenotypes(self):
        effect_sizes = self._get_effect_sizes()
        covs, edit_rates = _sample_all_from_data()

        guide_info["coverage"] = covs.astype(int)
        guide_info["edit_rate"] = edit_rates
        
        phenotypes = self.guide_info.reindex(self.guide_info.index.repeat(self.guide_info.coverage)).reset_index()
        phenotypes['edited'] = scipy.stats.bernoulli.rvs(phenotypes.edit_rate)
        phenotypes['reporter_edited'] = scipy.stats.bernoulli.rvs(phenotypes.edit_rate)
        phenotypes['phenotype_mean'] = phenotypes['edited']*phenotypes['effect_size']
        phenotypes['phenotype'] = scipy.stats.norm.rvs(loc = phenotypes.phenotype_mean, scale = 1)
        self.phenotypes = phenotypes



    def get_sorted_cell_counts(self):
        """Sort cells (phenotypes) into designated sorting scheme: bulk, 1/2/3/4"""
        if self.sorting_mode == "bins":
            self.get_sorted_cell_counts_quantiles(q = 5, bin_labels = self.bin_labels)
        elif self.sorting_mode == "topbot":
            self.get_sorted_cell_counts_quantiles(q = [0, 0.3, 0.7, 1], bin_labels = self.bin_labels)
        

    def get_sorted_cell_counts_quantiles(
        self, 
        q = [0, 0.2, 0.4, 0.6, 0.8, 1],
        bin_labels = None
        ):
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
                raise ValueError("provided length of bin_labels \
                    doesn't match the number of bins made by quantile")      

        phenotype_df["sorted_bin"] = pd.qcut(
            phenotype_df.phenotype,
            q = q,
            labels = bin_labels
        )

        guide_info = phenotype_df[["target_id", "guide_num"]].drop_duplicates()

        counts = _aggregate_counts(
            phenotype_df, 
           group_by = ["target_id", "guide_num", "sorted_bin"]
        ).reset_index()
        
        bulk = phenotype_df.sample(n = n_bulk_cells, replace = True)
        bulk_counts = _aggregate_counts(bulk, "bulk")

        sorted_counts = guide_info.merge(counts, on = ["target_id", "guide_num"])
        sorted_counts = sorted_counts.pivot(
            index = ["target_id", "guide_num"], 
            columns = "sorted_bin",
            values = ["guide", "target_edit", "reporter_edit"])

        sorted_counts.columns = ['_'.join(col[::-1]).strip() for col in sorted_counts.columns.values]
        sorted_counts = sorted_counts.reset_index().merge(
            bulk_counts, on = ["target_id", "guide_num"])
        sorted_counts = sorted_counts.fillna(0)
        return(sorted_counts)


    def get_sorted_read_counts(self, cell_counts : pd.DataFrame):
        measures = ["guide", "target_edit", "reporter_edit"]
        samples = self.sorting_bins + ["bulk"]
        read_counts = cell_counts[["target_id", "guide_num"]].reset_index(
            drop = True)
        for sample in samples:
            all_guide_count = cell_counts["{}_guide".format(sample)].sum()
            for measure in measures:
                bin_measure = "{}_{}".format(sample, measure)
                cell_counts["{}_p".format(sample)] = cell_counts[bin_measure] / all_guide_count
                assert not bin_measure in read_counts.columns
                read_counts[bin_measure] = \
                    scipy.stats.binom.rvs(
                        self.n_reads_per_sample, 
                        cell_counts["{}_p".format(sample)]
                    )
        return(read_counts)


    def simulate_rep(self):
        cell_counts = self.get_sorted_cell_counts()
        read_counts = self.get_sorted_read_counts(cell_counts)
        return(read_counts)


    def simulate_reps(self):       
        for rep in range(self.nreps):
            self.screen_res.append(simulate_rep())
            

