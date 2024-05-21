proc newRep { sel type color rep imol} {
    mol selection $sel
    mol representation $type
    mol addrep $imol
    mol showrep $imol $rep on
    mol modcolor $rep $imol $color
}

set dir "/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/08-Extract_molecules/MD_extract_g16_conformers"

display projection orthographic
axes location off
color Display Background white
display depthcue off

set listFiles {}
lappend listFiles [list $dir/MD_extract_018_gaussian_allign.mol2]
lappend listFiles [list $dir/MD_extract_011_gaussian_allign.mol2]
lappend listFiles [list $dir/MD_extract_000_gaussian_allign.mol2]
lappend listFiles [list $dir/MD_extract_001_gaussian_allign.mol2]
lappend listFiles [list $dir/MD_extract_005_gaussian_allign.mol2]
lappend listFiles [list $dir/MD_extract_013_gaussian_allign.mol2]
lappend listFiles [list $dir/MD_extract_014_gaussian_allign.mol2]


foreach ifile $listFiles {
    mol addfile $ifile type mol2
    set imol [molinfo top]
}
set imol_ref [molinfo top]
mol rename $imol_ref "OptimizedConformers"
mol delrep 0 $imol_ref
set rep1 0
newRep "all" "CPK" "Name" $rep1 $imol_ref
animate goto start

