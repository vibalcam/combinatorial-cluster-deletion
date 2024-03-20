# Combinatorial Approximations for Cluster Deletion: Simpler, Faster, and Better

Code for fast combinatorial methods to solve cluster deletion (DegMFP and RatMFP) and combinatorial solver for STC-LP.

Repo for the paper:
**Combinatorial Approximations for Cluster Deletion: Simpler, Faster, and Better**

## Code Organization

### Experiments

The experiments can be found in the `Experiments` folder.
`run-cd-functions.jl` contains the code to run a single instance of the cluster deletion problem with a given method.
The experiments are separated between `CD-experiments` and `merge-experiments`, with 

The `CD-experiments` folder contains the experiments for the cluster deletion problem.
- `run_many_main.jl` runs all the experiments from the paper
- `cd-table-print-results.jl` prints the results from the experiments in a table format and plots the figures
- `mpf-tests.jl` can be used to run the MFP methods faster

### Methods

You can find the implementations for the deterministic MFP methods (DegMFP and RatMFP) in `src/run_cd_functions.jl`.

The implementation for the combinatorial solver for STC-LP can be found in `src/cd_lp_minst_relaxations.jl`.

### Data

This repo contains some, but not all of the datasets used in the experiments. Other graphs may be accessed via the SNAP repository or the suitesparse matrix collection.

Details for how to download and standardize all snap graphs for experiments can be found in the folder data/snap-graphs, in file `standardize-snap-graphs.jl`.

## Acknowledgements

Our implementation is based on nveldt's repo [FastCC-via-STC](https://github.com/nveldt/FastCC-via-STC/tree/main).
