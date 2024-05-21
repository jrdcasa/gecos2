#!/bin/bash
#SBATCH --partition=all
#SBATCH --exclude=""
#SBATCH --cpus-per-task=4
#SBATCH --mem=2000M
#SBATCH --time=150:00:00
#SBATCH --job-name=C25H28O13_Conf_Cluster_0005

g16legacy_root=/optnfs/gaussian/g16_legacy
GAUSS_SCRDIR=$TMPDIR
source $g16legacy_root/bsd/g16.profile
export g16legacy_root GAUSS_SCRDIR
g16 C25H28O13_Conf_Cluster_0005.com
