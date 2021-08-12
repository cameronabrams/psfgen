# stripsplit.tcl (c) 2019-2021 cameron f abrams cfa22@drexel.edu
#
# Given a PSF file, and mutually exclusive atomselections "A" and "B", 
# strip all atoms not in A or B and write PSF and stripped DCDs for each
# of A, B and the union of A and B.
#

proc usage {} {
   puts "vmd -dispdev text -e stripsplit.tcl -args -Asel <str-w-dashes-instead-of-spaces>"
   puts "    -Bsel <str-w-dashes-instead-of-spaces> -psf <xxx.psf> -pdb <xxx.pdb> \\ "
   puts "    [-stride <int>] [-Aname <str-w-dashes-instead-of-spaces>] \\ "
   puts "    [-Bname <str-w-dashes-instead-of-spaces>] [-Aname <str-w-dashes-instead-of-spaces>] \\ "
   puts "    <dcd1> <dcd2> ..."
}

proc have_common_elements { list1 list2 } {
   foreach i $list1 j $list2 {
      if { $i in $list2 } { return 1 }
      if { $j in $list1 } { return 1 }
   }
   return 0
}

if {![info exists PSFGEN_BASEDIR]} {
  # see if user set an environment variable
  if {[info exists env(PSFGEN_BASEDIR)]} {
      set PSFGEN_BASEDIR $env(PSFGEN_BASEDIR)
  } else {
      set PSFGEN_BASEDIR $env(HOME)/research/psfgen
  }
}

source $PSFGEN_BASEDIR/scripts/vmdrc.tcl
package require pbctools

# define atomselect strings for A and B
set A_selstr ""
set B_selstr ""
set A_name "A"
set B_name "B"
set AB_name "AB"
set stride 1
set MKPSF 1
set PSF ""
set PDB ""
set DCDLIST [list]
# get the replica number and stride from the command-line
#puts "[llength $argv] $argv"
for { set i 0 } { $i < [llength $argv] } { incr i } {
   puts "$i [lindex $argv $i]"
    if { [lindex $argv $i] == "-stride" } {
       incr i
       set stride [lindex $argv $i]
    } elseif { [lindex $argv $i] == "-Asel" } {
       incr i
       set A_selstr [join [split [lindex $argv $i] "-"] " "]
    } elseif { [lindex $argv $i] == "-Bsel" } {
       incr i
       set B_selstr [join [split [lindex $argv $i] "-"] " "]
    } elseif { [lindex $argv $i] == "-Aname" } {
       incr i
       set A_name [join [split [lindex $argv $i] "-"] " "]
    } elseif { [lindex $argv $i] == "-Bname" } {
       incr i
       set B_name [join [split [lindex $argv $i] "-"] " "]
    } elseif { [lindex $argv $i] == "-ABname" } {
       incr i
       set AB_name [join [split [lindex $argv $i] "-"] " "]
    } elseif { [lindex $argv $i] == "-psf" } {
       incr i
       set PSF [lindex $argv $i]
#    } elseif { [lindex $argv $i] == "-pdb" } {
#       incr i
#       set PDB [lindex $argv $i]
    } elseif { [lindex $argv $i] == "-h" } {
       usage
       exit
    } else {
       lappend DCDLIST [lindex $argv $i]
    }
}
#puts "Aselstr $A_selstr Bselstr $B_selstr"
if { ! [file exists $PSF] } {
   puts "Error: psf $PSF not found."
   exit
}
#if { ! [file exists $PDB] } {
#   puts "Error: $PDB not found."
#   exit
#}
foreach dcd $DCDLIST {
   if { ! [file exists $dcd] } {
      puts "Error: dcd $dcd not found."
      exit
   }
}

package require psfgen
# load the overall simulation PSF file into psfgen
readpsf $PSF

# load the base molecule topology and coordinates
mol load psf $PSF ;#pdb $PDB
# determine list of unique segids 
set all_segids [lsort -unique [[atomselect top all] get segid]]

set Asel [atomselect top $A_selstr]
set Anum [$Asel num]
if { $Anum == 0 } {
   puts "Error: atomselection \"$A_selstr\" has no atoms."
   exit
}
set Bsel [atomselect top $B_selstr]
set Bnum [$Bsel num]
if { $Bnum == 0 } {
   puts "Error: atomselection \"$B_selstr\" has no atoms."
   exit
}

set Asegs [lsort -unique [$Asel get segid]]
set Bsegs [lsort -unique [$Bsel get segid]]
# assert: Asegs and Bsegs have no common elements
if { [have_common_elements $Asegs $Bsegs] } {
   puts "Error: A and B have segids in common."
   puts "I'm bailing out!"
   exit
}

set Xsegs [list]
foreach s $all_segids {
   if {[lsearch -exact $Asegs $s] >= 0} {
      # s is one of A's segids
      puts "found $s in Asegs"
      set inA [atomselect top "segid $s and ($A_selstr)"]
      set notinA [atomselect top "segid $s and not ($A_selstr)"]
      if {[$notinA num] == 0} {

      } else {
         foreach rn [$notinA get resid] {
            puts "appending $s $rn to Xsegs"
            lappend Xsegs [list $s $rn]
         }
      }
   } elseif {[lsearch -exact $Bsegs $s] >= 0} {
      # s is one of B's segids
      puts "found $s in Bsegs"
      set inB [atomselect top "segid $s and ($B_selstr)"]
      set notinB [atomselect top "segid $s and not ($B_selstr)"]
      if {[$notinB num] == 0} {

      } else {
         foreach rn [$notinB get resid] {
            puts "appending $s $rn to Xsegs"
            lappend Xsegs [list $s $rn]
         }
      }
   } else {
      puts "appending $s to Xsegs"
      lappend Xsegs [list $s]
   }
}

puts "Removing all atoms not in A or B: $Xsegs"
foreach xs $Xsegs {
   set tmpstr [join $xs " "]
   puts "Executing delatom $tmpstr"
   if {[llength $xs]==1} {
      delatom $xs
   } else {
      delatom [lindex $xs 0] [lindex $xs 1]
   }
}
# what remains is the PSF of the complex; write it
writepsf "${AB_name}.psf"

# delete the ligand segid; what remains is the PSF of the target; write it
foreach as $Asegs {
   delatom $as
}
writepsf "${B_name}.psf"
resetpsf

# read back in the complex psf; delete target segids; what remains is the ligand PSF; write it
readpsf "${AB_name}.psf"
foreach bs $Bsegs {
   delatom $bs
}
writepsf "${A_name}.psf"
resetpsf

mol delete top

# generate PDB files for each type of system from the
# input coordinates in first frame of first dcd
mol new $PSF
set working_id [molinfo top get id]
foreach dcd $DCDLIST {
   puts "Reading $dcd ..."
   animate read dcd $dcd skip $stride waitfor all
}
set ABsel [atomselect top "($A_selstr) or ($B_selstr)"]
set Asel [atomselect top "($A_selstr)"]
set Bsel [atomselect top "($B_selstr)"]
#pbc unwrap -all -sel "segid $targ_seg $lig_seg"
set aln_work [atomselect $working_id "($A_selstr) or ($B_selstr)"]
set aln_do [atomselect $working_id "($A_selstr) or ($B_selstr)"]
set aln_ref [atomselect $working_id "($A_selstr) or ($B_selstr)"]
$aln_ref frame 0
for { set i 0 } { $i < [molinfo top get numframes]} { incr i } {
    $aln_work frame $i
    $aln_do frame $i
    $aln_do move [measure fit $aln_work $aln_ref]
}
$ABsel frame 0
$ABsel writepdb "${AB_name}.pdb"
$Asel frame 0
$Asel writepdb "${A_name}.pdb"
$Bsel frame 0
$Bsel writepdb "${B_name}.pdb"
animate write dcd "${AB_name}.dcd" waitfor all sel $ABsel $working_id
animate write dcd "${A_name}.dcd" waitfor all sel $Asel $working_id
animate write dcd "${B_name}.dcd" waitfor all sel $Bsel $working_id

exit

