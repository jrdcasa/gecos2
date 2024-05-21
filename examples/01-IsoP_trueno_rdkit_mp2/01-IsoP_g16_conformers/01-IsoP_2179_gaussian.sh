#!/bin/bash
#SBATCH --partition=simacro
#SBATCH --exclude="trueno36, trueno37, trueno38, trueno59"
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000M
#SBATCH --job-name=01-IsoP_2179_gaussian

g16legacy_root=/opt/gaussian/g16_legacy
GAUSS_SCRDIR="$TMPDIR"
source $g16legacy_root/bsd/g16.profile
export g16legacy_root GAUSS_SCRDIR
$g16legacy_root/g16 01-IsoP_2179_gaussian.com
