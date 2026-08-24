#!/bin/bash
#SBATCH --job-name=eegglm
#SBATCH --output=/home/msequestro/logs/eegglm_%A_%a.log
#SBATCH --error=/home/msequestro/logs/eegglm_%A_%a.log
#SBATCH --partition=firstgen,secondgen,newgen,lastgen,dellgen
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-15:00:00
#SBATCH --array=1-50
#SBATCH --mail-user=matteo.sequestro@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --exclude=node43

module purge
unset LD_LIBRARY_PATH
unset SSL_CERT_FILE
unset SSL_CERT_DIR

module load MATLAB/R2024a

export MATLABROOT=$(dirname $(dirname $(which matlab)))
export LD_LIBRARY_PATH=$MATLABROOT/sys/os/glnxa64

RUN_ID=${SLURM_ARRAY_TASK_ID}

matlab -nodisplay -nosplash -singleCompThread -r \
  "tsglm_ssh_fit_eeg($RUN_ID); exit"