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

source $PSFGEN_BASEDIR/scripts/vmdrc.tcl

set PSF none.psf
set INPUTNAME none-input
set GMXTOPNAME none.top
set centerselstr "protein or glycan"
set PSFOUT tg_reordered.psf
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
  if { $arg == "-psf-out" } {
    incr a
    set PSFOUT [lindex $argv $a]
  }
  if { $arg == "--center-sel-str" } {
    incr a
    set centerselstr [lindex $argv $a]
  }
}

if { ! [file exists $PSF ] } {
   puts "Error: Input file $PSF not found."
   exit
}
if { ! [file exists ${INPUTNAME}.coor ] } {
   puts "Error: Input file ${INPUTNAME}.coor not found."
   exit
}
if { ! [file exists ${INPUTNAME}.xsc ] } {
   puts "Error: Input file ${INPUTNAME}.xsc not found."
   exit
}
if { [string compare $GMXTOPNAME "none.top"] == 0 } {
   puts "Error: You must specify the name of the output gromacs topology using -top <name>.top"
   exit
}
if { ! [file exists $parinp ] } {
   puts "Error: Input file $parinp not found."
   exit
}

# Read the par.inp generated by do_py.sh.  This file is a set of NAMD configuration statements like:
# parameters <parameter-file-name>
# We take the second column (1) as the names of local CHARMM parameter files to feed to writegmxtop
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

# We have to make sure connected fragments have identical chain IDs and segnames.
# However, when cfaparsepdb.py builds using psfgen, it gives unique segnames to
# each glycan.  This, among other things, must be undone.
# Note: If the parent PDB file gives a glycan a fully unique chain ID, cfaparsepdb.py
# currently (8/9/21) by default does not revert that to the chain ID of the protein to 
# which it is connected.  So, those "orphan" glycans need to be processed separately.
set protein_chains [lsort -unique [[atomselect top "protein and name CA"] get chain]]
puts "# Protein chains detected: $protein_chains"
set protein_segnames [lsort -unique [[atomselect top "protein and name CA"] get segname]]
puts "# Protein segnames detected: $protein_segnames"
foreach P $protein_chains {
    set oions [atomselect top "chain $P and ion"]
    set oions_count [$oions num]
    if { $oions_count > 0 } {
       $oions set chain I
       $oions set segname ION
       puts "# Changing chain to I and segname to ION for $oions_count ions in chain $P"
    }
    set glycans [atomselect top "chain $P and glycan"]
    set glycans_count [$glycans num]
    if { $glycans_count > 0 } {
      $glycans set segname $P
      puts "# Changing glycan segname to $P"
    }
    set xw [atomselect top "chain $P and water"]
    set xw_count [$xw num]
    if { $xw_count > 0 } {
      $xw set chain W
      $xw set segname WTX
      puts "# Changing crystal water chain to W and segname to WTX"
    }
}
# Here we identify and process orphan chains.  First, we identify all chainID's that have not 
# already been identified as having protein atoms; these are the orphans.  Then we 
# interrogate each orphan chain to find out what *existing* protein chain that exactly
# one of its atoms is bonded to.  If such a connection exists, we then set the orphan chain's ID 
# and segname to the chain ID of that "outside" protein atom.
set non_protein_chains [lsort -unique [[atomselect top "not protein and not ion and not water"] get chain]]
puts "# Chains with non-protein/ion/water atoms detected: $non_protein_chains"
set non_protein_segnames [lsort -unique [[atomselect top "not protein and not ion and not water"] get segname]]
puts "# Segnames with non-protein/ion/water atoms detected: $non_protein_segnames"
set orphan_chains [list]
foreach npc $non_protein_chains {
  if { $npc in $protein_chains } {
  } else {
    lappend orphan_chains $npc
  }
}
if { [llength $orphan_chains] > 0 } {
  puts "# Chains that have not been processed for topogromacs: $orphan_chains"
  foreach OC $orphan_chains {
    set ocsel [atomselect top "chain $OC"]
    set ocatoms [$ocsel get index]
    set ocbl [$ocsel getbonds]
    set outside [list]
    foreach abl $ocbl {
      foreach abli $abl {
        if { $abli in $ocatoms } {

        } else {
          lappend outside $abli
        }
      }
    }
    if { [llength $outside] > 0 } {
      puts "# Orphan chain $OC: atom(s) outside of chain to which atom(s) in chain are bonded: $outside"
      set ownerchain [lsort -unique [[atomselect top "index [lindex $outside 0]"] get chain]]
      if { $ownerchain in $protein_chains } {
        puts "#   -> Setting chain ID and segname to $ownerchain"
        $ocsel set chain $ownerchain
        $ocsel set segname $ownerchain 
      }
    } else {
      puts "# Orphan chain $OC is independent and not connected to any protein atom."
    }
  }
}

puts "# Reanalyzing topology after necessary chain/segname changes"
mol reanalyze top
puts "# Done reanalyzing"

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
animate write psf $PSFOUT $newmol
animate write pdb tg_reordered.pdb $newmol

mol delete 0
mol delete $newmol

mol new $PSFOUT
mol addfile tg_reordered.pdb

set all [atomselect top all]
pbc set $cell

set fp [open "tg-cell-nm.in" "w"]
puts $fp "[expr [lindex $cell 0] / 10.0]  [expr [lindex $cell 1] / 10.0] [expr [lindex $cell 2] / 10.0]"
close $fp

topo writegmxtop $GMXTOPNAME $CHARMMPARFILES 
pbc wrap -now -compound fragment -centersel $centerselstr -center com
[atomselect top "all"] writepdb tg_needsbox.pdb

puts "# Generated final outputs $PSFOUT and $GMXTOPNAME"
puts "# Intermediate files tg_needsbox.pdb and tg-cell-nm.in need further processing through gmx editconf."

exit
