import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils import *



def load_and_filter_data(sim_type, introduce_after, sim_length, parameters_set, h_threshold):
    """
    Charge et filtre les données pour un type de simulation donné.

    Parameters:
    sim_type (str): Le type de simulation (ex. 'ssd', 'ra', 'sss_densest').
    introduce_after (int): Temps après lequel commencer l'analyse.
    sim_length (int): Durée totale de la simulation.
    parameters_set (int): Le set de paramètres utilisé pour la simulation.
    h_threshold (float): Le seuil d'hétéroplasmie.

    Returns:
    pd.DataFrame: DataFrame filtré avec les données de la simulation.
    """
    h_df = pd.read_pickle(f"parameters_set{parameters_set}_{sim_type}_h_{100*h_threshold:.0f}.pkl")
    h_df = h_df[h_df["cell"] == -1]  # Consider only the aggregated data across all cells
    filtered_h_df = h_df[(h_df["t"] >= introduce_after) & (h_df["t"] <= sim_length)]
    return filtered_h_df

if __name__ == "__main__":
    # Parameters
    introduce_after = 2500  # Modifier selon le moment où vous souhaitez commencer
    sim_length = 5000    # Modifier selon la durée de la simulation souhaitée
    parameters_set = 6
    h_threshold = 1.0

    # Simulations types
    sim_types = ["ra", "sss_densest","ssd"]
    colors = ['blue', 'red', "purple"]
    
    #sim_types=["sss_densest"]
    #colors=['red']
    
    fig, ax = plt.subplots(figsize=(6,4))
    
    for sim_type, color in zip(sim_types, colors):
        filtered_h_df = load_and_filter_data(sim_type, introduce_after, sim_length, parameters_set, h_threshold)
        
        
        # Plot the increase in mutant fraction against time
        #ax.plot(filtered_h_df["t"], filtered_h_df["h_mean"], linewidth=0.8, color=color, label=sim_type)
        #ax.fill_between(filtered_h_df["t"], filtered_h_df["h_mean"] - filtered_h_df["h_sem"], filtered_h_df["h_mean"] + filtered_h_df["h_sem"], color=color, alpha=0.25)
        # Renommer 'sss_densest' en 'sss' pour la légende
        legend_label = "sss" if sim_type == "sss_densest" else sim_type
        
        # Tracer l'augmentation de la fraction mutante en fonction du temps
        ax.plot(filtered_h_df["t"], filtered_h_df["h_mean"], linewidth=0.8, color=color, label=legend_label)
        ax.fill_between(filtered_h_df["t"], filtered_h_df["h_mean"] - filtered_h_df["h_sem"], filtered_h_df["h_mean"] + filtered_h_df["h_sem"], color=color, alpha=0.25)
    
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\langle f \rangle$")
    ax.legend()
    plt.tight_layout()
    plt.show()
