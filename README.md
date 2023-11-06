# Pesto manuscript scripts

This is a code supplement intended to make the analyses and figures presented in the Pesto methodology manuscript reproducible for the readers and the reviewers.

## Simulation study

The simulation analyses are done in the following files:
* `scripts/simulate_trees.jl`: simulating trees under the birth-death-shift process
* `scripts/simulation_inference.jl`: performing the inference on the simulated trees. There are 3 inference settings, and 2000 simulated trees, so 6000 runs, meaning that it will be slow on a personal computer. There is a HPC script `scripts/simulation_inference_job.sh` with an example of how we ran our analyses on our cluster with 250 GB memory and 80 CPU cores in parallel. This script will write the results as extended newick strings (one to each file), as well as some summary matrices (JLD2 format) to the `simulated_trees` directory.
* `scripts/confusion_matrix.R`: reads all the completed analyses in R, and compares the true simulation history vs the inferred rate estimates and shift estimates. This calculates metrics such as the false positive rate, false negative rate, proportional error in branch-rate estimates, etc.

## Reproduce figures

Once these are ran, one can reproduce Figures 2-11 in the main text using the scripts prefixed 2-11. Figure 1 is a TiKz figure, and we don't add the code to actually plot it, but we include a script to calculate the probability densities along the branches of the three-taxon tree.

## How to run scripts

In order to run the scripts, you will need:
* An up-to-date version of R, with several packages (i.e. treeio, ggtree, tidytree, dplyr etc.)
* An up-to-date version of Julia, with the modules `Pesto` and `BirthDeathSimulation` installed, and some other modules for importing/exporting files, and making figures (`StatsPlots`and `CairoMakie`).
* Scripts are ran from the top-level of the repository, in the terminal. For example, to make Fig. 2, you would write: `julia scripts/02_quantiles.jl`

## Exact simulation reproducibility

The simulated trees will not be exactly the same as the ones shown in the manuscript. In order to get the exact same trees, you can download a `.zip` from this URL: [insert link, 500 MB]

