# topogromacs script

package require topotools
package require pbctools
mol new my_6m0j_i.psf; # created by psfgen
set inputname sol-stage2.restart
mol addfile ${inputname}.coor
pbc readxst ${inputname}.xsc

set cell [molinfo top get {a b c}]

# make some changes to undo my original segnames and chains selections
# 1. put the Zn ion into the ION segment
set zn [atomselect top "chain A and name ZN"]
$zn set chain I
$zn set segname ION

# 2. put the gycans into their parent chains' segmments
set as [atomselect top "chain A and segname AS"]
$as set segname A
set es [atomselect top "chain E and segname ES"]
$es set segname E

# 3. put the crystal waters into chain W and give them
# a common segname
set awx [atomselect top "chain A and segname AWX"]
$awx set chain W
$awx set segname WTX
set ewx [atomselect top "chain E and segname EWX"]
$ewx set chain W
$ewx set segname WTX

mol reanalyze top

# at this point fragments are out of order, so let's try to fix that
set fragsellist [list]
set bigsel [atomselect top "not (water or ions)"]
foreach frag [lsort -u [$bigsel get fragment]] {
  set fsel [atomselect top "fragment $frag"]
  lappend fragsellist $fsel
}
set othersel [atomselect top "(water or ions)"]
lappend fragsellist $othersel
set newmol [::TopoTools::selections2mol $fragsellist]
animate write psf my_6m0j_i_reordered.psf $newmol
animate write pdb my_6m0j_i_reordered.pdb $newmol

mol delete 0
mol delete $newmol

mol new my_6m0j_i_reordered.psf
mol addfile my_6m0j_i_reordered.pdb

set all [atomselect top all]
pbc set $cell

set fp [open "tg-cell-nm.in" "w"]
puts $fp "[expr [lindex $cell 0] / 10.0]  [expr [lindex $cell 1] / 10.0] [expr [lindex $cell 2] / 10.0]"
close $fp

#[atomselect top all] writepdb "tg_my_6m0j_i_step1.pdb"
topo writegmxtop tg_my_6m0j_i.top [list par_all36_carb.prm  par_all36_prot.prm toppar_all36_carb_glycopeptide.str  toppar_water_ions.str]
pbc wrap -now -compound fragment -centersel "protein or glycan" -center com
[atomselect top "all"] writepdb tg_my_6m0j_i_needsbox.pdb
puts "Generated tg_my_6m0j_i.top, tg_my_6m0j_i_needsbox.pdb, and tg-cell-nm.in"
exit
