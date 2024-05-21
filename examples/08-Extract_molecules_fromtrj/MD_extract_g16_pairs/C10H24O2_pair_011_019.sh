#!/bin/bash
#SBATCH --partition=all
#SBATCH --exclude=""
#SBATCH --cpus-per-task=8
#SBATCH --mem=8000M
#SBATCH --time=150:00:00
#SBATCH --job-name=C10H24O2_pair_011_019

hostname
g16root=/optnfs/gaussian
source $g16root/g16_legacy/bsd/g16.profile
g16 C10H24O2_pair_011_019.com
