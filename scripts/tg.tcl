# topogromacs script
# cameron f abrams cfa22@drexel.edu
if {![info exists PSFGEN_BASEDIR]} {
  # see if user set an environment variable
  if {[info exists env(PSFGEN_BASEDIR)]} {
      set PSFGEN_BASEDIR $env(PSFGEN_BASEDIR)
  } else {
      set PSFGEN_BASEDIR $env(HOME)/research/psfgen
  }
}

set PSF none.psf
set INPUTNAME none-input
set GMXTOPNAME none.top
set centerselstr "protein or glycan"
set CHARMMPARFILES [list]
set parinp par.inp
for { set a 0 } { $a < [llength $argv] } { incr a } {
  set arg [lindex $argv $a]
  if { $arg == "-psf" } {
    incr a
    set PSF [lindex $argv $a]
  }
  if { $arg == "-i" } {
    incr a
    set INPUTNAME [lindex $argv $a]
  }
  if { $arg == "-top" } {
    incr a
    set GMXTOPNAME [lindex $argv $a]
  }
  if { $arg == "-parinp" } {
    incr a
    set parinp [lindex $argv $a]
  }
  if { $arg == "--center-sel-str" } {
    incr a
    set centerselstr [lindex $argv $a]
  }
}

if { ! [file exists $PSF ] } {
   puts "Error: $PSF not found."
   exit
}
if { ! [file exists ${INPUTNAME}.coor ] } {
   puts "Error: ${INPUTNAME}.coor not found."
   exit
}
if { ! [file exists ${INPUTNAME}.xsc ] } {
   puts "Error: ${INPUTNAME}.xsc not found."
   exit
}
if { [string compare $GMXTOPNAME "none.top"] == 0 } {
   puts "Error: You must specify the name of the gromacs topology to create with -top <name>.top"
   exit
}
if { ! [file exists $parinp ] } {
   puts "Error: $parinp not found."
   exit
}

set pf [open $parinp "r"]
while {[gets $pf line] >= 0} {
  set tokens [split $line]
  lappend CHARMMPARFILES [lindex $tokens 1]
}
close $pf

package require topotools
package require pbctools
puts "# Reading $PSF, ${INPUTNAME}.coor, and ${INPUTNAME}.xsc"
mol new $PSF
mol addfile ${INPUTNAME}.coor
pbc readxst ${INPUTNAME}.xsc

set cell [molinfo top get {a b c}]

set protein_chains [lsort -unique [[atomselect top "protein and name CA"] get chain]]
set protein_segnames [lsort -unique [[atomselect top "protein and name CA"] get segname]]
foreach P $protein_chains {
    set oions [atomselect top "chain $P and ion"]
    $oions set chain I
    $oions set segname ION
    set glycans [atomselect top "chain $P and glycan"]
    $glycans set segname $P
    set xw [atomselect top "chain $P and water"]
    $xw set chain W
    $xw set segname WTX
}


# make some changes to undo my original segnames and chains selections
# 1. put the Zn ion into the ION segment
#set zn [atomselect top "chain A and name ZN"]
#$zn set chain I
#$zn set segname ION

# 2. put the gycans into their parent chains' segmments
#set as [atomselect top "chain A and segname AS AG00 AG01 AG02 AG03"]
#$as set segname A
#set es [atomselect top "chain E and segname ES EG00"]
#$es set segname E

# 3. put the crystal waters into chain W and give them
# a common segname
#set awx [atomselect top "chain A and segname AWX AW"]
#$awx set chain W
#$awx set segname WTX
#set ewx [atomselect top "chain E and segname EWX EW"]
#$ewx set chain W
#$ewx set segname WTX

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
animate write psf tg_reordered.psf $newmol
animate write pdb tg_reordered.pdb $newmol

mol delete 0
mol delete $newmol

mol new tg_reordered.psf
mol addfile tg_reordered.pdb

set all [atomselect top all]
pbc set $cell

set fp [open "tg-cell-nm.in" "w"]
puts $fp "[expr [lindex $cell 0] / 10.0]  [expr [lindex $cell 1] / 10.0] [expr [lindex $cell 2] / 10.0]"
close $fp

#[atomselect top all] writepdb "tg_my_6m0j_i_step1.pdb"
topo writegmxtop $GMXTOPNAME $CHARMMPARFILES 
#[list par_all36_carb.prm  par_all36_prot.prm toppar_all36_carb_glycopeptide.str  toppar_water_ions.str]
pbc wrap -now -compound fragment -centersel $centerselstr -center com
[atomselect top "all"] writepdb tg_needsbox.pdb
puts "Generated $GMXTOPNAME, tg_needsbox.pdb, and tg-cell-nm.in"
exit
