# Le Spectre de Fréquence des Sites Distingue les Mécanismes Derrière l'Accumulation des Mutations de l'ADNmt


## Structure du Dépôt
Le dépôt est organisé en plusieurs répertoires pour séparer les différentes composantes du projet, facilitant ainsi la navigation et la maintenance.

- `include/` : Contient les fichiers d'en-tête C qui déclarent les fonctions, macros, types de données et valeurs de paramètres utilisés dans les fichiers sources C.
- `src/` : Contient les fichiers sources C qui simulent la réplication et la dégradation des molécules d'ADNmt dans une cellule.
- `data_analysis/` : Contient les scripts Python utilisés pour analyser les données générées par les simulations.




### Exemple d'Utilisation
Si on veut effectuer pour quelques simulations à la main, mais le plus gros a tét fait avec les codes en .pbs
1. **Compiler les programmes C**
    - Naviguer dans le répertoire `src/`
    ```bash
    cd src
    ```
    - Compiler les fichiers sources
    ```bash
    gcc ra_single_sim.c lib_sim.c -o ra_single_sim.exe -lgsl -lgslcblas -lm
    ```
2. **Exécuter les simulations**
    ```bash
    ./ra_single_sim.exe 1
    ./ra_single_sim.exe 2
    ```
Ici, deux simulations avec les graines 1 et 2 sous le modèle de l'avantage réplicatif (RA) sont effectuées en utilisant les valeurs de paramètres spécifiées dans `include/parameters.h`. Un nouveau répertoire `src/parameters_set1` est créé, contenant les résultats des simulations `ra_sim_populations1.txt`, `ra_sim_site_frequency_spectrum1.txt`, `ra_sim_populations2.txt`, `ra_sim_site_frequency_spectrum2.txt`, ainsi qu'une liste des valeurs de paramètres `parameters.txt`.

### Analyse des Données
1. **Copier les données des simulations dans le répertoire pour l'analyse des données**
    ```bash
    cd
    cp src/parameters_set1/ra_sim*.txt data_analysis/
    cd data_analysis
    ```
2. **Compiler les données des simulations**
    ```bash
    python compile_data.py 1 ra pop
    python compile_data.py 1 ra sfs
    ```
    Les dataframes Pandas `parameters_set1_ra_pop_master.pkl` et `parameters_set1_ra_sfs_master.pkl` sont générés et sauvegardés. Ces dataframes contiennent les données compilées des 2 simulations.
3. **Combiner les données de populations et le spectre de fréquence des sites**
    ```bash
    python combine_pop_sfs.py 1 ra
    ```
    Le dataframe Pandas `parameters_set1_ra_combined_master.pkl` est généré et sauvegardé.
4. **Ensemble des données**
    ```bash
    python ensemble_data.py 1 ra 1.0
    ```
    Les dataframes Pandas `parameters_set1_ra_h_100.pkl` et `parameters_set1_ra_sfs_100.pkl` contenant la fraction de mutants et le spectre de fréquence des sites à travers les simulations sont générés. Ces dataframes sont généralement de petite taille et peuvent être copiés sur votre machine locale.
5. **Visualiser les données**
    - Tracer la fraction de mutants en fonction du temps
    ```bash
    python example_plot_h.py
    ```
    - Tracer le spectre de fréquence des sites
    ```bash
    python example_plot_sfs.py
    ```
6. **Effectuer des tests statistiques**

    Supposons qu'en plus des simulations sous le modèle RA, les étapes ci-dessus soient également suivies pour simuler, compiler, et rassembler les données sous le modèle de survie stochastique du plus dense (SSD), obtenant ainsi le dataframe Pandas `parameters_set1_ssd_sfs_100.pkl`. Ensuite, la commande suivante peut être exécutée pour effectuer un test de Mann-Whitney U sur la différence des spectres de fréquence des sites sous les modèles RA et SSD.
    ```bash
    python example_mann-whitney.py
    ```

