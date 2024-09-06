#!/usr/bin/env sh
#SBATCH --job-name=bears_SCM
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=b.kopperud@lmu.de
#SBATCH --mem=2GB
#SBATCH --output=logs/BDS-scm-%a.log
#SBATCH --error=logs/BDS-scm-%a.err
#SBATCH --qos=normal_prio
#SBATCH --ntasks=1

module load gnu openmpi

#source "scripts/env.sh"

REPLICATE=$(cat scripts/arg_list.txt | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')

#mpirun -np 2 rb-mpi --args ${TIMESLICES} --file scripts/bears_scm.Rev
rb --args ${REPLICATE} --file scripts/primates_scm.Rev
