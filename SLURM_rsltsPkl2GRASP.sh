#!/usr/local/bin/bash
#SBATCH --job-name=ABI_0
#SBATCH --nodes=1
#SBATCH --time=0:15:00
#SBATCH -o log/output.%A-%a
#SBATCH -e log/error.%A-%a
#SBATCH --array=0
#SBATCH --account=s2417

date
hostname
echo "---Running Simulation N="${SLURM_ARRAY_TASK_ID}"---"
python convert_rsltsPkl_to_groupedGraspRun.py
wait
exit 0
