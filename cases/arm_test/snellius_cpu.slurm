#!/bin/bash
#SBATCH --job-name=arm

# Log file `stdout` and `stderr` (%j is replaced with a unique job ID):
#SBATCH --output=mhh-%j.out
#SBATCH --error=mhh-%j.err

# New (2023) Genoa CPU partition (https://servicedesk.surfsara.nl/wiki/display/WIKI/Snellius+usage+and+accounting):
#SBATCH --partition=genoa

# Genoa nodes have 192 cores/node. Slurm automatically determines the required number of nodes:
#SBATCH -n 64
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-core=1

# Maximum wall clock time:
#SBATCH -t 00:59:00

# E-mail notifications, options are: {none, begin, end, fail, all}:
#SBATCH --mail-type=end
#SBATCH --mail-user=swar@wetsus.nl

module purge
module load 2022
module load CMake/3.23.1-GCCcore-11.3.0
module load foss/2022a
module load netCDF/4.9.0-gompi-2022a
module load CUDA/11.8.0
module load Clang/13.0.1-GCCcore-11.3.0

export OMPI_MCA_fcoll="two_phase"
export OMPI_MCA_io_ompio_bytes_per_agg="512MB"

srun ./microhh init arm
srun ./microhh run arm