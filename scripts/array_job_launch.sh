#-- Launch slurm job array
echo -e "Launching analyses... "
sbatch -p krypton --array="1-$(wc -l scripts/arg_list.txt | cut -d ' ' -f1)" ./scripts/job_scm.sh
echo "done"
