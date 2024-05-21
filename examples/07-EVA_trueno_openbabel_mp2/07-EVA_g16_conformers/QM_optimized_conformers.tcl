proc newRep { sel type color rep imol} {
    mol selection $sel
    mol representation $type
    mol addrep $imol
    mol showrep $imol $rep on
    mol modcolor $rep $imol $color
}

set dir "/home/jramos/PycharmProjects/GeCos/examples/07-EVA_trueno_openbabel_mp2/07-EVA_g16_conformers"

display projection orthographic
axes location off
color Display Background white
display depthcue off

set listFiles {}
lappend listFiles [list $dir/07-EVA_001_gaussian_allign.mol2]
lappend listFiles [list $dir/07-EVA_021_gaussian_allign.mol2]


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

