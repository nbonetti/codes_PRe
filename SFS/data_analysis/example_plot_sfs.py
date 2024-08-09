import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils import *

if __name__ == "__main__":
    # Parameters
    parameters_set = 6
    h_threshold = 1.0
    #sim_types=["ra"]
    #colors=['red']
    sim_types = ["ra",  "sss_densest","ssd"]
    colors = ['blue', 'red',"purple"]

    fig, ax = plt.subplots(figsize=(6,4))
    
    for sim_type, color in zip(sim_types, colors):
        # Load data
        sfs_df = pd.read_pickle(f"parameters_set{parameters_set}_{sim_type}_sfs_{100*h_threshold:.0f}.pkl")
        sfs = sfs_df[sfs_df["cell"] == -1]["sfs"][0]  # Consider only the aggregated data across all cells
        sfs = remove_homoplasmy(sfs)
        
        # Plot the site frequency spectrum
        n_bins = 20  # set number of bins
        plot_sfs(ax, sfs, n_bins, density=True, label=sim_type, type="line", marker="o", markersize=3, color=color)
    
    ax.set_xlabel("h")
    ax.set_ylabel("density")
    ax.legend()
    plt.tight_layout()
    plt.show()



