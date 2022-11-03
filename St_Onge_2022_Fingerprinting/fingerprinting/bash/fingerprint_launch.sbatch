#!/bin/bash

#SBATCH --job-name=fp_launch
#SBATCH --output=%x-%a.out                             # Write standard output and standard error to this file [uses job name (%x) and job ID (%A)]
#SBATCH --array=1-7
#SBATCH --time=10:00:00                                # Time limit in hrs:min:sec before job is killed (strongly affects how quickly job is started)
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2                              # Number of CPU cores per task; 1 unless you know you can use multiple cores
#SBATCH --mem=4G                                       # Memory per node
#SBATCH --account=CC_USER
#SBATCH --mail-user=EMAIL    # Specify your email address for notifications
#SBATCH --mail-type=ALL                                # Email you about these events (options: NONE, BEGIN, END, FAIL, ALL)
# you can comment out any unnecessary SBATCH options with an additional #

# print working directory, hostname, date and Job ID
pwd; hostname; date; echo "Job ID: ${SLURM_JOB_ID}"; echo "Task ID: ${SLURM_TASK_ID}"

# Setting up the variables for the command
OUTPUT_DIR="/ABSOLUTE/PATH/TO/OUTPUT/output_folder" 
PYT_SCRIPT="/ABSOLUTE/PATH/TO/SCRIPT/fingerprinting.py"
PATH_FILE="/ABSOLUTE/PATH/TO/FILE/CONTAINING/ARGUMENTS/file.txt"
ARGS_PERMUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${PATH_FILE})

#Note: Path file refers to a text file containing the different permutations of arguments
# we want the fingerprint to run. It needs to be built by the user, and it modifies
# the number of "slurm_arrays" needed to input above.

#Loading Python and virtual environment
module load python/3.8
source ~/python_envs/fp_project/bin/activate 

#Preparing command to launch script
PYT_CMD="python -u ${PYT_SCRIPT} -t within -p=Schaefer -c=pcor --output=${OUTPUT_DIR} ${ARGS_PERMUT}"

# Print the Python command for reference
echo ${PYT_CMD}

# Launch eval of the Python command.
eval ${PYT_CMD}