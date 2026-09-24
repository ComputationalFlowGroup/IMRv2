#!/bin/bash

# Copyright (C) Mauro Rodriguez
# SPDX-License-Identifier: GPL-3.0-or-later

# Request an hour of runtime:
#SBATCH --time=6:00:00

# Controls the minimum/maximum number of nodes allocated to the job
#SBATCH -N 1
#SBATCH -p batch
#SBATCH -c 32
# Default resources are 1 core with 2.8GB of memory.

# Use more memory (4GB):
#SBATCH --mem=60G

# Specify a job name:
#SBATCH -J forsim

# Specify an output file
#SBATCH -o forsim-%j.out
#SBATCH -e forsim-%j.err

# Run f_multirun, in the same directory as this batch script, and save its results.
matlab-threaded -nodisplay -r "[newtonian, neohook, voigt, ucm, oldb] = f_multirun(); save('multirun.mat'); exit"

