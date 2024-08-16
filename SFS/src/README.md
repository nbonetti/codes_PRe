# SRC

This directory contains the C source files to simulate mtDNA mutations under different models.

## Files

- **_debug.c**  
  Contains debugging and printing functions. It can be compiled with other source files to create an executable, which can then be run with `gdb` for debugging.

- **lib_sim.c**  
  Provides functions for the general simulation framework.
    - Implements a structure called `Dict`, similar to the Python `dict` type, to store site frequency spectrum information.
    - Contains functions to update the system's state based on replication, degradation, or diffusion events. The occurrence of these events is probabilistically determined according to the system's propensity. The rule for updating wildtype propensity (`wildtype_propensity_update`) is provided, while the mutants' propensity is determined in the `*single_sim.c` file under the specified model.
    - Also contains functions to write data to text files.

- **lib_sim_2.c**  
  Similar to `lib_sim.c`, but specific to the `sss_densest` model, where mutants are only slowed down (`delta=1`).

- **ra_single_sim.c**  
  Source for simulating a system under the RA model.

- **ssd_single_sim.c**  
  Source for simulating a system under the SSD model.

- **sss_densest_single_sim.c**  
  Source for simulating a system under the SSS model, with `delta=1` as in our models, though this value can be modified if dealing with the densest+slowest mutants case.

## Usage
GCC is required. Without debugging, compile the `*single_sim.c` source file with `lib_sim.c`:
```bash
gcc [ra, ssd, sss_densest, or wildtype]_single_sim.c lib_sim.c -o single_sim.exe -lgsl -lgslcblas -lm
./single_sim.exe 1
