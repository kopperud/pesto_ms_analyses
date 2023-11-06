#!/usr/bin/env sh
#SBATCH --job-name=bdshift_power_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b.kopperud@lmu.de
#SBATCH --mem=250GB
#SBATCH --output=logs/inference_bdshift.log
#SBATCH --error=logs/inference_bdshift.err
#SBATCH --qos=normal_prio_large
#SBATCH --ntasks=80
#SBATCH --nodes=1
#SBATCH --partition=lemmium

module load R/4.2.3 gnu openblas

export R_HOME="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R"
export LD_LIBRARY_PATH="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R/lib"
#julia scripts/setup/dependencies.jl > logs/jldep.txt
#Rscript scripts/setup/dependencies.R > logs/Rdep.txt
julia --threads 80 scripts/simulation_inference.jl > output/screen.txt
