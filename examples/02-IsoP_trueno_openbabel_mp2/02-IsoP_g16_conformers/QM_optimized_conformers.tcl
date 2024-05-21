proc newRep { sel type color rep imol} {
    mol selection $sel
    mol representation $type
    mol addrep $imol
    mol showrep $imol $rep on
    mol modcolor $rep $imol $color
}

set dir "/home/jramos/PycharmProjects/GeCos/examples/02-IsoP_trueno_openbabel_mp2/02-IsoP_g16_conformers"

display projection orthographic
axes location off
color Display Background white
display depthcue off

set listFiles {}
lappend listFiles [list $dir/02-IsoP_028_gaussian_allign.mol2]
lappend listFiles [list $dir/02-IsoP_036_gaussian_allign.mol2]
lappend listFiles [list $dir/02-IsoP_007_gaussian_allign.mol2]
lappend listFiles [list $dir/02-IsoP_012_gaussian_allign.mol2]
lappend listFiles [list $dir/02-IsoP_014_gaussian_allign.mol2]
lappend listFiles [list $dir/02-IsoP_013_gaussian_allign.mol2]
lappend listFiles [list $dir/02-IsoP_000_gaussian_allign.mol2]
lappend listFiles [list $dir/02-IsoP_001_gaussian_allign.mol2]
lappend listFiles [list $dir/02-IsoP_030_gaussian_allign.mol2]
lappend listFiles [list $dir/02-IsoP_008_gaussian_allign.mol2]


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

