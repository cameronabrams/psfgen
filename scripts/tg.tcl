# topogromacs vmd script
# cameron f abrams cfa22@drexel.edu
# 2021

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
set PDBOUT tg_intermediate.pdb
set ALTNAMDCONFIG ""
set CELLDIMFILE "tg-cell-nm.in"
set CHARMMPARFILES [list]
set parinp par.inp
set default_resid_pad 100
for { set a 0 } { $a < [llength $argv] } { incr a } {
  set arg [lindex $argv $a]
  if { $arg == "-psf" } {
    incr a
    set PSF [lindex $argv $a]
  } elseif { $arg == "-i" } {
    incr a
    set INPUTNAME [lindex $argv $a]
  } elseif { $arg == "-top" } {
    incr a
    set GMXTOPNAME [lindex $argv $a]
  } elseif { $arg == "-parinp" } {
    incr a
    set parinp [lindex $argv $a]
  } elseif { $arg == "-opsf" } {
    incr a
    set PSFOUT [lindex $argv $a]
  } elseif { $arg == "-opdb" } {
    incr a
    set PDBOUT [lindex $argv $a]
  } elseif { $arg == "--alt-namd-config" } {
    incr a
    set ALTNAMDCONFIG [lindex $argv $a]
  } elseif { $arg == "--cell-dim-file" } {
    incr a
    set CELLDIMFILE [lindex $argv $a]
  } elseif { $arg == "--center-sel-str" } {
    incr a
    set centerselstr [lindex $argv $a]
  } else {
    puts "Warning: option $arg not recognized; ignoring."
  }
}

if { ! [file exists $PSF ] } {
   vmdcon -err "Input file $PSF not found."
   exit
}
if { ! [file exists ${INPUTNAME}.coor ] } {
   vmdcon -err "Input file ${INPUTNAME}.coor not found."
   exit
}
if { ! [file exists ${INPUTNAME}.xsc ] } {
   vmdcon -err "Input file ${INPUTNAME}.xsc not found."
   exit
}
if { [string compare $GMXTOPNAME "none.top"] == 0 } {
   vmdcon -err "You must specify the name of the output gromacs topology using -top <name>.top"
   exit
}
if { ! [file exists $parinp ] } {
   vmdcon -info "Input file $parinp containing NAMD \"parameters xxx\" config directives not found."
   if { ! [file exists $ALTNAMDCONFIG] } {
     vmdcon -err "Could not find a NAMD config file to extract parameter names from"
     exit
   } else {
     vmdcon -info "Generating $parinp from $ALTNAMDCONFIG"
     exec grep "^parameters $ALTNAMDCONFIG" > $parinp
   }
}

# Read the par.inp generated by do_py.sh.  This file is a set of NAMD configuration statements like:
# parameters <parameter-file-name>.  If a namd configuration exists, one can generate by redirecting
# output of "grep ^parameters my_file.namd" to "par.inp".  We take the second column (1) as the names 
# of local CHARMM parameter files to feed to writegmxtop
set pf [open $parinp "r"]
while {[gets $pf line] >= 0} {
  set tokens [split $line]
  lappend CHARMMPARFILES [lindex $tokens 1]
}
close $pf

package require topotools
package require pbctools
vmdcon -info "Reading $PSF, ${INPUTNAME}.coor, and ${INPUTNAME}.xsc"
mol new $PSF
mol addfile ${INPUTNAME}.coor
pbc readxst ${INPUTNAME}.xsc

set cell [molinfo top get {a b c}]

# We have to make sure connected fragments have identical chain IDs and segnames.
# However, when cfaparsepdb.py builds using psfgen, it gives unique segnames to
# each glycan.  This, among other things, must be undone.
# Note: If the parent PDB file gives a glycan a fully unique chain ID, cfaparsepdb.py
# So, those "orphan" glycans need to be processed separately.
set protein_chains [lsort -unique [[atomselect top "protein and name CA"] get chain]]
vmdcon -info "Protein chains detected: $protein_chains"
set protein_segnames [lsort -unique [[atomselect top "protein and name CA"] get segname]]
vmdcon -info "Protein segnames detected: $protein_segnames"
set all_ions [atomselect top "ion"]
set ion_chains [list]
set ion_segnames [list]
set fixI 0
if { [$all_ions num] > 0 } {
  set ion_chains [lsort -unique [[atomselect top "ion"] get chain]]
  vmdcon -info "Ion chains detected: $ion_chains"
  foreach ic $ion_chains {
    if { $ic in $protein_chains } {
      vmdcon -info "Shared protein/ion chain ID: $ic"
      if { $ic == "I" } { 
        vmdcon -warn "Warning: your topology has a protein chain \"I\"."
        vmdcon -warn "Chain \"I\" is reserved for initially free ions."
        vmdcon -warn "The protein chain \"I\" will be assigned a new chainID."
        set fixI 1
      }
    }
  }
  set ion_segnames [lsort -unique [$all_ions get segname]]
  vmdcon -info "Ion segnames detected: $ion_segnames"
  if { [llength $ion_segnames] > 1} {
    vmdcon -info "All ions will be given segname \"ION\""
  }
  $all_ions set segname "ION"
}
foreach P $protein_chains {
  # handle case that there are ions sharing a chain id with a protein chain;
  # just change segname to ION, leave chain ID alone
  set oions [atomselect top "chain $P and ion"]
  set oions_count [$oions num]
  if { $oions_count > 0 } {
    $oions set segname ION
    vmdcon -info "Changing segname to \"ION\" for $oions_count ions in protein-chain $P"
  }
  # handle glycans sharing chainID with a protein; set their segname to that of the protein
  set glycans [atomselect top "chain $P and glycan"]
  set glycans_count [$glycans num]
  if { $glycans_count > 0 } {
    $glycans set segname $P
    vmdcon -info "Changing glycan segname to $P"
  }
  # handle any crystal waters; just put them all in the same chain with unique segname
  set xw [atomselect top "chain $P and water"]
  set xw_count [$xw num]
  if { $xw_count > 0 } {
    $xw set chain W
    $xw set segname WX${P}
    vmdcon -info "Changing crystal water chain to W and segname to WX${P}"
  }
}
# Here we identify and process orphan chains.  First, we identify all chainID's that have not 
# already been identified as having protein atoms; these are the orphans.  Then we 
# interrogate each orphan chain to find out what *existing* protein chain that exactly
# one of its atoms is bonded to.  If such a connection exists, we then set the orphan chain's ID 
# and segname to the chain ID of that "outside" protein atom.
set non_protein_chains [lsort -unique [[atomselect top "not protein and not ion and not water"] get chain]]
vmdcon -info "Chains with non-protein/ion/water atoms detected: $non_protein_chains"
set non_protein_segnames [lsort -unique [[atomselect top "not protein and not ion and not water"] get segname]]
vmdcon -info "Segnames with non-protein/ion/water atoms detected: $non_protein_segnames"
set orphan_chains [list]
foreach npc $non_protein_chains {
  if { $npc in $protein_chains } {
  } else {
    lappend orphan_chains $npc
  }
}
if { [llength $orphan_chains] > 0 } {
  vmdcon -info "Chains that have not been processed for topogromacs: $orphan_chains"
  foreach OC $orphan_chains {
    set ocsel [atomselect top "(not protein and not ion and not water) and chain $OC"]
    set ocatoms [$ocsel get index]
    set ocbl [$ocsel getbonds]
    set outside [list]
    foreach abl $ocbl {
      foreach abli $abl {
        if { $abli in $ocatoms } {
          # do nothing
        } else {
          lappend outside $abli
        }
      }
    }
    if { [llength $outside] > 0 } {
      vmdcon -info "Orphan chain $OC: atom(s) outside of chain to which atom(s) in chain are bonded: $outside"
      set ownerchain [lsort -unique [[atomselect top "index [lindex $outside 0]"] get chain]]
      if { $ownerchain in $protein_chains } {
        set highest_resid [lindex [lsort -u -integer -decreasing [[atomselect top "chain $ownerchain and not water"] get resid]] 0]
        set this_offset [expr $highest_resid + $default_resid_pad]
        $ocsel set chain $ownerchain
        $ocsel set segname $ownerchain
        set bad_resids [$ocsel get resid]
        set new_resids [list]
        foreach bi $bad_resids {
          lappend new_resids [expr $bi + $this_offset]
        }
        $ocsel set resid $new_resids
        vmdcon -info "  -> Setting chain ID and segname to $ownerchain with lowest resid [lindex $new_resids 0]"
      }
    } else {
      vmdcon -info "Orphan chain $OC is independent and not connected to any protein atom."
    }
  }
}

# put ions BACK into chain I if they aren't there -- this can happen if there is an orphan/protein chain I already
# rename any protein/glycan chain I to the next available chain ID
set pci [atomselect top "(protein or glycan) and chain I"]
set pcinum [$pci num]
if { $pcinum > 0 } {
  # need a new chain ID for this chain
  set all_chains [lsort -unique [[atomselect top all] get chain]]
  set available_chainIDs [list]
  foreach ch [split "ABCDEFGHIJKLMNOPQRSTUVWXYZ" {}] {
    if { $ch in $protein_chains } {
      vmdcon -info "searching for new chain ID for old chain I... $ch is taken."
    } else {
      vmdcon -info "searching for new chain ID for old chain I... $ch is available."
      lappend available_chainIDs $ch
    }
  }
  vmdcon -info "Available chain IDs: [join $available_chainIDs ',']"
  set newch [lindex $available_chainIDs 0]
  vmdcon -info "selected $newch"
  $pci set chain $newch
  $pci set segname $newch
} else {
  vmdcon -info "No protein/glycan atoms with chainID I found."
}
set ii [atomselect top "ion"]
if { [$ii num] > 0 } {
  set iicids [lsort -unique [$ii get chain]]
  vmdcon -info "Ions chainID(s) detected: [join $iicids ',']"
  set inoti [atomselect top "ion and not chain I"]
  if { [$inoti num] > 0 } {
    vmdcon -info "Setting chain ID of all atoms selected by macro \"ion\" to I"
    $inoti set chain "I"
  }
  # Paranoid sanity check -- macro "ion" and "segname ION" should select identical atoms
  # but ions may have chain IDs other than I if they are associated with a protein chain
  set is [atomselect top "segname ION"]
  set iindx [$ii get index]
  set isndx [$is get index]
  vmdcon -info "ion macro indices: $iindx"
  vmdcon -info "segname ION indices: $isndx"
  set diff [lmap n [concat $iindx $isndx] {
    # Skip the elements that are in both lists
    if {$n in $iindx && $n in $isndx} continue
    set n
  }]
  if { [llength $diff] > 0 } {
    vmdcon -info "Error: some ions selected by macro ion differ from ions selected using segname ION"
    vmdcon -info "This means the code needs fixing."
    vmdcon -info "ion macro indices: $iindx"
    vmdcon -info "segname ION indices: $isndx"
    exit
  }
}

vmdcon -info "Reanalyzing topology after necessary chain/segname changes"
mol reanalyze top
vmdcon -info "Done reanalyzing"

# at this point fragments are out of order, so let's try to fix that
vmdcon -info "Using TopoTools to reorder out-of-order fragments"
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

# generate cell dimensions in nanometers
set fp [open $CELLDIMFILE "w"]
puts $fp "[expr [lindex $cell 0] / 10.0]  [expr [lindex $cell 1] / 10.0] [expr [lindex $cell 2] / 10.0]"
close $fp

vmdcon -info "Calling writegmxtop..."
vmdcon -info "Command: topo writegmxtop $GMXTOPNAME $CHARMMPARFILES"
topo writegmxtop $GMXTOPNAME $CHARMMPARFILES 
vmdcon -info "Done with call to writegmxtop.  Centering and writing PDB."
pbc wrap -now -compound fragment -centersel $centerselstr -center com
[atomselect top "all"] writepdb $PDBOUT

vmdcon -info "Generated final output files $PSFOUT and $GMXTOPNAME"
vmdcon -info "Next command:"
vmdcon -info "> gmx editconf -f $PDBOUT -o tg_final.pdb -box `cat $CELLDIMFILE`"
exit
