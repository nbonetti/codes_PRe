# SRC

Ce répertoire contient les fichiers sources en C pour simuler les mutations de l'ADNmt sous différents modèles.

## Fichiers

- **_debug.c**  
  Contient des fonctions de débogage et d'impression. Peut être compilé avec d'autres fichiers sources pour créer un exécutable, que vous pouvez ensuite exécuter avec `gdb` pour le débogage.

- **lib_sim.c**  
  Fournit des fonctions pour le cadre général de simulation.
    - Implémente une structure appelée `Dict`, semblable au type dict en Python, pour stocker les informations du spectre de fréquence des sites.
    - Contient des fonctions pour mettre à jour l'état du système en fonction des événements de réplication, dégradation ou diffusion. La survenue de ces événements est déterminée probabilistiquement selon la propension de l'état. La règle de mise à jour de la propension des types sauvages `wildtype_propensity_update` est fournie, tandis que celle des mutants est déterminée dans le fichier `*single_sim.c` sous le modèle spécifié.
    - Contient également des fonctions pour écrire les données dans des fichiers texte.
 
- **lib_sim_2.c**
  Même principe que pour celui d'avant mais celui-ci est spécifique au modèle sss_densest où les mutants sont uniquements ralentits 'delta=1'.

- **ra_single_sim.c**  
  Source pour simuler un système sous le modèle RA.

- **ssd_single_sim.c**  
  Source pour simuler un système sous le modèle SSD.

- **sss_densest_single_sim.c**  
  Source pour simuler un système sous le modèle SSS.


## Utilisation
GCC est requis. Sans débogage, compilez le fichier source `*single_sim.c` avec `lib_sim.c` :
```bash
gcc [ra, ssd, sss_densest, ou wildtype]_single_sim.c lib_sim.c -o single_sim.exe -lgsl -lgslcblas -lm
./single_sim.exe 1

