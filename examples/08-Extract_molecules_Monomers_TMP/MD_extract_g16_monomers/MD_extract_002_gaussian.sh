#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --exclude=""
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000M
#SBATCH --job-name=MD_extract_002_gaussian

g16legacy_root=/opt/g16
GAUSS_SCRDIR="$TMPDIR"
source $g16legacy_root/bsd/g16.profile
export g16legacy_root GAUSS_SCRDIR
$g16legacy_root/g16 MD_extract_002_gaussian.com
