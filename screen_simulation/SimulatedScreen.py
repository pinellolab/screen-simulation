from typing import List, Optional, Union
from inspect import ismethod

import numpy as np
import pandas as pd
from pprint import pprint
import scipy.stats


edit_rates_saved_data = np.loadtxt("/data/pinello/PROJECTS/2021_08_ANBE/data/07+1021_ANBE_data/crisprep_count/bulk_edit_rates.csv")
guide_coverage_saved_data = np.loadtxt("/data/pinello/PROJECTS/2021_08_ANBE/data/07+1021_ANBE_data/crisprep_count/bulk_coverage_sums.csv")


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
        n_targets: int = 700,
        n_guides_per_target: int = 5,
        n_total_cells: int = 10**6*8,
        n_bulk_cells = None,
        sorting_mode: str = "topbot",
        nreps:int = 4,
        n_genome_per_sample:int = int(10e4),
        n_reads_per_sample:int = 10**6*2,
        effect_size_distribution: Optional[Union[np.ndarray, str]] = None,
        guide_coverage_distribution: Optional[Union[np.ndarray, str]] = "data",
        edit_rate_distribution: Optional[Union[np.ndarray, str]] = "data",
        has_reporter = False,
        ):

        self.sorting_mode = sorting_mode
        if self.sorting_mode == "topbot":
            self.sorting_bins = ["bot", "mid", "top"]
        elif self.sorting_mode == "bins":
            self.sorting_bins =  ["bot", "low", "mid", "high", "top"]
        else:
            raise ValueError("Invalid sorting_mode: {}".format(sorting_mode))
        self.sample_names = self.sorting_bins + ["bulk"]
        self.sorting_bins_tested = [sb for sb in self.sorting_bins if sb != "mid"]

        self.n_reps = nreps
        self.n_targets = n_targets
        self.n_guides_per_target = n_guides_per_target
        self.n_guides = self.n_targets * self.n_guides_per_target

        if edit_rate_distribution is None:
            self.edit_rate_dist = np.array([1])
        elif edit_rate_distribution == "data":
            self.edit_rate_dist = edit_rates_saved_data
        elif isinstance(edit_rate_distribution, np.ndarray):
            self.edit_rate_dist = edit_rate_distribution
        else:
            raise ValueError("Invalid edit_rate_distribution: {}".format(edit_rate_distribution))

        if guide_coverage_distribution is None:
            self.guide_coverage_dist = np.array([1])
        elif guide_coverage_distribution == "data":
            self.guide_coverage_dist = guide_coverage_saved_data
        elif isinstance(guide_coverage_distribution, np.ndarray):
            self.guide_coverage_dist = guide_coverage_distribution
        else:
            raise ValueError("Invalid guide_coverage_distribution: {}".format(guide_coverage_distribution))

        if effect_size_distribution is None:
            self.effect_size_dist = np.concatenate((np.repeat(0, 60), np.linspace(0.1, 1, 10)))
        elif isinstance(effect_size_distribution, np.ndarray):
            self.effect_size_dist = effect_size_distribution
        else:
            raise ValueError("Invalid effect_size_distribution: {}".format(effect_size_distribution))

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


    def __repr__(self):
        out_str = "SimulatedScreen object with \n"
        for k, v in getattr(self, "__dict__").items():
            if k.startswith("_"): continue
            if ismethod(v): continue
            if isinstance(v, pd.DataFrame): continue
            if k in ["samples", "screen_res"]: continue
            out_str += "\t{} = {}\n".format(k, v)
        return(out_str)

    def get_id(self):
        params = "{}_{}_{}_{}".format(
            self.sorting_mode, 
            self.n_cells_per_rep, 
            self.n_reps, 
            self.n_reads_per_sample
        )
        return(params)


    def _sample_all_from_data(self):
        # samples (guide_count, edit_rate) from self.count_edit_stats
        edit_rate = np.random.choice(self.edit_rate_dist, size =self.n_guides, replace = True)
        guide_count_raw = np.random.choice(self.guide_coverage_dist, size = self.n_guides, replace = True)
        guide_count_norm = np.floor(guide_count_raw / self.guide_coverage_dist.sum() * self.n_total_cells)
        return((guide_count_norm, edit_rate))

    def _get_effect_sizes(self):
        effect_sizes = np.resize(self.effect_size_dist, self.n_targets)
    
        self.guide_info = pd.DataFrame({
            "target_id" : np.repeat(list(range(self.n_targets)), self.n_guides_per_target), 
            "guide_num" : np.tile(list(range(self.n_guides_per_target)), self.n_targets), 
            "effect_size" : np.repeat(effect_sizes, self.n_guides_per_target)
            })

        self.target_info = self.guide_info[["target_id", "effect_size"]].drop_duplicates().reset_index()
        return(effect_sizes)


    def simulate_cells(self):
        """ Simulates cells as dataframe as self.cells. """

        # Simulate the effect sizes of the variants.
        effect_sizes = self._get_effect_sizes()
        # Get coverage and edit rate pairs for each guide
        covs, edit_rates = self._sample_all_from_data()

        # Cell coverage and edit rate per guide
        self.guide_info["coverage"] = covs.astype(int)
        self.guide_info["edit_rate"] = edit_rates
        
        cells = self.guide_info.reindex(self.guide_info.index.repeat(self.guide_info.coverage)).reset_index()
        cells['edited'] = scipy.stats.bernoulli.rvs(cells.edit_rate)
        cells['phenotype_mean'] = cells['edited']*cells['effect_size']
        cells['phenotype'] = scipy.stats.norm.rvs(loc = cells.phenotype_mean, scale = 1) 
        #TODO: for scRNA readout this line should be different @Lucas

        if self.has_reporter:
            cells['reporter_edited'] = scipy.stats.bernoulli.rvs(cells.edit_rate)

        self.cells = cells

    def expand_cells(self):
        self.cells = NotImplemented #TODO

    def select_cells(self):
        self.cells = NotImplemented #TODO

    def sort_cells(self):
        """Sort cells into designated sorting scheme: bulk, 1/2/3/4"""
        if self.sorting_mode == "bins":
            self.sort_cells_quantiles(q = [0, 0.2, 0.4, 0.6, 0.8, 1], bin_labels = self.sorting_bins)
        elif self.sorting_mode == "topbot":
            self.sort_cells_quantiles(q = [0, 0.3, 0.7, 1], bin_labels = self.sorting_bins)
        

    def sort_cells_quantiles(
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

        self.cells["sorted_bin"] = pd.qcut(
            self.cells.phenotype,
            q = q,
            labels = bin_labels
        )



    def get_genomes(self):
        """
        Samples the genomes for each sequencing samples prepped.
        Number of genomes sampled is sampled as the multinomial 
        distribution with p proportional to cell counts per guide.

        Assigns self.samples attribute which is Dict[str -> pd.DataFrame]
        that stores sampled genomes per sample with UMI
        """
        self.samples = {} # Dict[str -> pd.DataFrame]
        cell_record_to_keep = ["target_id", "guide_num", "edited"]
        if self.has_reporter: cell_record_to_keep.append("reporter_edited")
        
        # sorted bins

        for sorting_bin, df in self.cells.groupby("sorted_bin", as_index = False):
            if self.n_genome_per_sample > len(df): 
                # All genomes are sampled
                self.samples[sorting_bin] = df[cell_record_to_keep]
            else:
                self.samples[sorting_bin] = df[cell_record_to_keep].sample(
                    n = self.n_genome_per_sample, replace = False)

        # bulk sample
        self.samples["bulk"] = self.cells[cell_record_to_keep].sample(
            n = self.n_genome_per_sample, replace = False)

        # Assign UMI
        for sample_name in self.sample_names:
            self.samples[sample_name]["UMI"] = list(range(len(self.samples[sample_name])))

        
    def amplify_reads(self):
        for sample_name, df in self.samples.items():
            amplify_count = np.ones(len(df)) #TODO: different amplification per molecule (UMI)
            amplify_count *= self.n_reads_per_sample / amplify_count.sum()
            self.samples[sample_name] = pd.DataFrame(np.repeat(df.values, amplify_count.astype(int), axis = 0))
            self.samples[sample_name].columns = df.columns

    def get_read_counts(self):
        agg_fn = {"edited":["count", "sum"]}
        if self.has_reporter:
            agg_fn["reporter_edited"] = "sum"
        measures_umi = list(map(lambda s: s + "_UMI", self.measures))
        samples_counts = self.guide_info
        for sample_name, df in self.samples.items():
            counts = df.groupby(["target_id", "guide_num"]).agg(agg_fn)
            counts_umi = df.drop_duplicates().groupby(["target_id", "guide_num"]).agg(agg_fn)
            counts.columns = ["{}_{}".format(sample_name, s) for s in self.measures]
            counts_umi.columns = ["{}_{}".format(sample_name, s) for s in measures_umi]
            samples_counts = samples_counts.merge(counts, on = ["target_id", "guide_num"]).merge(
                counts_umi, on = ["target_id", "guide_num"])
        return(samples_counts)

    def simulate_rep(self):
        self.simulate_cells()
        # TODO: growing cells and adapting for viability screen
        # self.expand_cells()
        # self.select_cells()
        self.sort_cells()
        self.get_genomes()
        self.amplify_reads()
        read_counts = self.get_read_counts()
        return(read_counts)


    def simulate_reps(self):       
        for rep in range(self.n_reps):
            self.screen_res.append(self.simulate_rep())
            

