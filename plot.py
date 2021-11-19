import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

def plot_compare_schemes(sens_df: pd.DataFrame, measure, cov, nrep, nreads, ax, title = ""):
    df_plt = sens_df.loc[(sens_df["measure"] == measure) &
                         (sens_df["cov"] == cov) &
                         (sens_df["nrep"] == nrep) &
                         (sens_df["nreads"] == nreads) , :]
    ymax = min(df_plt[["top", "high", "low", "bot"]].max().max() + 0.1, 1)
    df_plt_bins = df_plt.loc[(df_plt.sorting_mode == "bins"), :]
    df_plt_topbot = df_plt.loc[(df_plt.sorting_mode == "topbot"), :]

    ax.plot(df_plt_bins.effect_size, df_plt_bins.top, label = "4bins:top", c = "red")
    ax.plot(df_plt_bins.effect_size, df_plt_bins.high, label = "4bins:high", c = "orange")
    ax.plot(df_plt_bins.effect_size, df_plt_bins.low, label = "4bins:low", c = "skyblue")
    ax.plot(df_plt_bins.effect_size, df_plt_bins.bot, label = "4bins:bot", c = "blue")
    ax.plot(df_plt_bins.effect_size, df_plt_bins["any"], label = "4bins:any", c = "black")

    ax.plot(df_plt_topbot.effect_size, df_plt_topbot.top, label = "2bins:top", c = "red", linestyle = '--')
    ax.plot(df_plt_topbot.effect_size, df_plt_topbot.bot, label = "2bins:bot", c = "blue", linestyle ='--')
    ax.plot(df_plt_topbot.effect_size, df_plt_topbot["any"], label = "2bins:any", c = "black", linestyle ='--')
    
    
    ax.set_ylim((-0.05, 1.05))
    ax.set_xlabel("effect_size")
    ax.set_ylabel("sensitivity")
    ax.set_title(title)
    ax.legend()
    
def plot_compare_schemes_readnorm(sens_df: pd.DataFrame, measure, cov, nrep, ax, title = ""):
    """Temporary fn to plot nread = 5*10**5 for "bins" and 10**6 for "topbot"
    """
    nreads_bins = 5*10**5
    nreads_topbot = 10**6
    df_plt = sens_df.loc[(sens_df["measure"] == measure) &
                         (sens_df["cov"] == cov) &
                         (sens_df["nrep"] == nrep) &
                         (sens_df["nreads"].map(lambda x: x in [nreads_bins, nreads_topbot])) , :]
    ymax = min(df_plt[["top", "high", "low", "bot", "any"]].max(
    numeric_only = True).max() + 0.1, 1.05)
    df_plt_bins = df_plt.loc[(df_plt.sorting_mode == "bins")&
                         (df_plt["nreads"] == nreads_bins) , :]
    df_plt_topbot = df_plt.loc[(df_plt.sorting_mode == "topbot")&
                         (df_plt["nreads"] == nreads_topbot) , :]

    ax.plot(df_plt_bins.effect_size, df_plt_bins.top, label = "4bins:top", c = "red")
    ax.plot(df_plt_bins.effect_size, df_plt_bins.high, label = "4bins:high", c = "orange")
    ax.plot(df_plt_bins.effect_size, df_plt_bins.low, label = "4bins:low", c = "skyblue")
    ax.plot(df_plt_bins.effect_size, df_plt_bins.bot, label = "4bins:bot", c = "blue")
    ax.plot(df_plt_bins.effect_size, df_plt_bins["any"], label = "4bins:any", c = "black")

    ax.plot(df_plt_topbot.effect_size, df_plt_topbot.top, label = "2bins:top", c = "red", linestyle = '--')
    ax.plot(df_plt_topbot.effect_size, df_plt_topbot.bot, label = "2bins:bot", c = "blue", linestyle ='--')
    ax.plot(df_plt_topbot.effect_size, df_plt_topbot["any"], label = "2bins:any", c = "black", linestyle ='--')
    
    ax.set_ylim((-0.05, 1.05))
    ax.set_xlabel("effect_size")
    ax.set_ylabel("sensitivity")
    ax.set_title(title)
    ax.legend()
    
def plot_compare_measures(sens_df: pd.DataFrame, sorting_mode, cov, nrep, nreads, ax, title = ""):
    df_plt = sens_df.loc[(sens_df["sorting_mode"] == sorting_mode) &
                         (sens_df["cov"] == cov) &
                         (sens_df["nrep"] == nrep) &
                         (sens_df["nreads"] == nreads) , :]
    if sorting_mode == "bins":
        sorting_bins = ["top", "high", "low", "bot", "any"]
        cols = ["red", "orange", "skyblue", "blue", "black"]
    else:
        sorting_bins = ["top", "bot", "any"]
        cols = ["red", "blue", "black"]
    ymax = min(df_plt[sorting_bins].max().max() + 0.1, 1.05)
    df_plt_guide = df_plt.loc[(df_plt.measure == "guide"), :]
    df_plt_target = df_plt.loc[(df_plt.measure == "target_edit"), :]
    df_plt_reporter = df_plt.loc[(df_plt.measure == "reporter_edit"), :]

    for i, sorting_bin in enumerate(sorting_bins):
        ax.plot(df_plt_guide.effect_size, df_plt_guide[sorting_bin], 
                label = "guide_{}".format(sorting_bin), c = cols[i])
        ax.plot(df_plt_target.effect_size-0.01, df_plt_target[sorting_bin],
               label = "target_{}".format(sorting_bin), c = cols[i], 
               linestyle = "dashed")
        ax.plot(df_plt_reporter.effect_size-0.02, df_plt_reporter[sorting_bin],
               label = "reporter_{}".format(sorting_bin), c = cols[i], 
               linestyle = "dotted")

    ax.set_ylim((-0.05, 1.05))
    ax.set_xlabel("effect_size")
    ax.set_ylabel("sensitivity")
    ax.set_title(title)
    #ax.legend()
    