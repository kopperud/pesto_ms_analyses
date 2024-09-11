#!/usr/bin/env sh
#SBATCH --job-name=bds_sim_inference
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b.kopperud@lmu.de
#SBATCH --mem=800GB
#SBATCH --output=logs/inference_bdshift.log
#SBATCH --error=logs/inference_bdshift.err
#SBATCH --qos=low_prio_res
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=200
#SBATCH --partition=krypton

module load R/4.2.3 gnu openblas

export R_HOME="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R"
export LD_LIBRARY_PATH="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R/lib"
#julia scripts/setup/dependencies.jl > logs/jldep.txt
#Rscript scripts/setup/dependencies.R > logs/Rdep.txt
echo ${SLURM_CPUS_PER_TASK} > logs/ntasks.txt

julia --threads ${SLURM_CPUS_PER_TASK} scripts/simulation_inference.jl > output/screen.txt
