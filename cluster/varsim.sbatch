#!/bin/bash
#SBATCH --job-name postVar # Set a name for your job. This is useful if you have multiple jobs queued.
#SBATCH --partition short    # Slurm partition to use
#SBATCH --ntasks 1         # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 0-8:00       # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=1024    # Memory limit for each tasks (in MB)
#SBATCH -o varsim.out     # File to which STDOUT will be written
#SBATCH -e varsim.out
#SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=amitmeir@uw.edu # Email to which notifications will be sent
#SBATCH --array=2-99

module load R

R CMD BATCH --no-save --no-restore "--args seed='${SLURM_ARRAY_TASK_ID}'" varMSsim_A.R