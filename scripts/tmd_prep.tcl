# 
# The purpose of this script is to direct VMD to generate a pdb file
# for a targeted MD simulation that starts with a system defined using
# a PSF/PDB or PSF/COOR and targets a system defined by a PDB.
#
# The resulting target PDB is congruent with (same number of atoms in same order as)
# the input PSF/PDB/COOR pair.  
#
# The two basic steps are
# 1. align the target PDB to the system PDB using some alignment basis
# 2. change coordinates of a selection of atoms in the system using 
#    coordinates in the target PDB
#
# In addition to identifying the three input files (the system PSF and PDB or COOR,
# and the target PDB), the user must specify the alignment basis selection and the 
# target selection.  
#
# ALIGNMENT BASIS.  This is a comma-separated list of colon-delimited pairs.  
# Each pair indicates unique sets of congruent residues in the system and 
# target.  A unique set of residues is expressed in the form X_xxx-yyy, where
# "X" is the chainID, xxx is the N-terminal-most residue in the range and 
# yyy is the C-terminal-most residue in the range.  Each pair is checked
# to ensure the two selections are indeed congruent up to backbone atoms.
# The list of pairs is used to generate a single atomselection string that when
# used in either the system molecule or the target molecule will produce
# a congruent set of residues permitting a [measure fit] command.
#
# TARGET SELECTION.  This is also a comma-separated list of colon-delimited
# pairs.  After alignment of the target onto the system via the 
# alignment basis, the system atoms in the target selection have their
# coordinates overwritten by the corresponding coordinates in the target.
#
# A full selection of all atoms in the system is then output as the
# TMD PDB.
#
# Cameron Abrams cfa22@drexel.edu

# generates an atomselect string from the command-line shortcode
# shortcode format is ABC_xxx-yyy/zzz-ggg/....
# where the characters in the substring to the left of the "_" delimiter
# are chain IDs, and the substring after the "_" delimiter 
# is a slash-delimited list (no spaces!) of one or more residue ranges.
# Each residue-range is a dash-delimited list of resid's. If a dash is 
# omitted, the single resid is inferred.
proc make_atsel { str } {
    set fp [split $str "_"]
    set cids [join [split [lindex $fp 0] {}] " "]
    set rrl [split [lindex $fp 1] "/"]
    set nrr [llength $rrl]
    set rrs "(chain $cids and (resid"
    for { set i 0 } { $i < $nrr } { incr i } {
        set lr [split [lindex $rrl $i] "-"]
        if { [llength $lr] == 1 } {
            set rrs "$rrs [lindex $lr 0]"
        } else {
            set rrs "$rrs [lindex $lr 0] to [lindex $lr 1]"
        }
    }
    set rrs "$rrs))"
    return $rrs
}

proc build_atom_select_strings { csl } {
    set pairs [split $csl ","]
    set ll [list]
    set rl [list]
    set np [llength $pairs]
    for { set i 0 } { $i < $np } { incr i } {
        set pair [split [lindex $pairs $i] ":"]
        set ld [lindex $pair 0]
        set rd [lindex $pair 1]
        lappend ll [make_atsel $ld]
        lappend rl [make_atsel $rd]
    }
    set return_me [list [join $ll " or "] [join $rl " or "]]
    return $return_me
}

proc mapresids { shortcode &map } {
    upvar 1 ${&map} arr
    set myl [split $shortcode " "]
    for {set i 0} {$i<[llength $myl]} {incr i} {
	    set rawwrd [lindex $myl $i]
	    set newwrd $rawwrd
	    set chrs [split $rawwrd ""]
        set anum ""
	    set fst -1
	    set lst 0
	    for {set j 0} {$j<[llength $chrs]} {incr j} {
	        set c [lindex $chrs $j]
	        if {[string is integer $c]} {
                set anum "$anum$c"
		        if {$fst==-1} {
		            set fst $j
		        }
		        set lst $j
            }
        }
	    if { $anum != "" } {
	        set newwrd [string replace $rawwrd $fst $lst $arr($anum)]
	        lset myl $i $newwrd
        }
    }
    set newshortcode [join $myl " "]
    return $newshortcode
}

proc main { argv } {
    set seland "protein and backbone"
    set argc [llength $argv]
    for { set i 0 } { $i < $argc } { incr i } {
        if { [lindex $argv $i] == "--alignment-basis" } {
            incr i
            set albs [lindex $argv $i]
        }
        if { [lindex $argv $i] == "--target-sel" } {
            incr i
            set tass [lindex $argv $i]
        }
        if { [lindex $argv $i] == "--system" } {
            incr i
            set res [split [lindex $argv $i] ","]
            set psf [lindex $res 0]
            set cor [lindex $res 1]
        } 
        if { [lindex $argv $i] == "--target-pdb" } {
            incr i
            set tarpdb [lindex $argv $i]
        }
        if { [lindex $argv $i] == "--output-pdb" } {
            incr i
            set outpdb [lindex $argv $i]
        }
        if { [lindex $argv $i] == "--seland" } {
            incr i
            set seland [lindex $argv $i]
        }
        if { [lindex $argv $i] == "--map21" } {
            incr i
	        # argument is name of tcl script that defines
	        # an associative array called "mymap"
            source [lindex $argv $i]
	        if { ! [info exists mymap] } {
                puts "Error: [lindex $argv $i] does not define mymap.  No target resid renumbering is done."
	        }
        }
    }
    
    mol new $psf
    mol addfile $cor
    set sys [molinfo top get id]
    mol new $psf
    mol addfile $cor
    set savsys [molinfo top get id]
    # zero occupancies for all but the target atoms
    set a1 [atomselect $sys "protein and backbone"]
    set a2 [atomselect $savsys "protein and backbone"]
    puts "[$a1 num] [$a2 num]"
    [atomselect $sys all] set occupancy 0
    $a1 set occupancy 1

    mol new $tarpdb 
    set tar [molinfo top get id]
    set alb [build_atom_select_strings $albs]
    set tas [build_atom_select_strings $tass]

    # Resid mapping:  in general it is assumed that the resid's in the
    # shortcodes are in the same numbering scheme; in other words, the
    # system and target have the same numbering.  If they do not,
    # it should be assumed that the resids in the system shortcodes (which
    # are alwasy to the LEFT of the colon) follow system numbering and
    # resids in the target shortcodes follow the target numbering.  It is
    # up the user to be sure the atomselections implied by the two
    # shortcodes are congruent.
    #
    # In some cases in which the system and target have different numbering,
    # and the user wishes to express all shortcodes in one numbering, the user
    # must also supply a "mapping" so that atomselections built from the
    # shortcodes indicate the correct resids for their respective molecules.
    # For example, the shortcode "A_101-110:A_101-110" indicates resids 101 to 110
    # on chain A both the system and target.  However, if this numbering 
    # applies to the target but not to the system, the resids 101-110 on chain A
    # will be translated to their corresponding numberings in the system by 
    # the user-supplied map.
    if { [info exists mymap ] } {
	    set mapped [mapresids [lindex $alb 0] mymap]
	    lset alb 0 $mapped
	    set mapped [mapresids [lindex $tas 0] mymap]
	    lset tas 0 $mapped
    }
    puts "System: $psf $cor"
    puts "Target pdb: $tarpdb"
    puts "Alignment basis: SYS([lindex $alb 0]) and TARG([lindex $alb 1])"
    puts "Target selection: SYS([lindex $tas 0]) and TARG([lindex $tas 1])"

    # perform an alignment-based move of the entire target
    set tarall [atomselect $tar all]
    set taralb [atomselect $tar "$seland and ([lindex $alb 1])"]
    set sysalb [atomselect $sys "$seland and ([lindex $alb 0])"]
    $tarall move [measure fit $taralb $sysalb]

    # copy target coordinates onto the system
    set sysrcv [atomselect $sys "$seland and ([lindex $tas 0])"]
    set syssavx [$sysrcv get x]
    set syssavy [$sysrcv get y]
    set syssavz [$sysrcv get z]
    set targiv [atomselect $tar "$seland and ([lindex $tas 1])"]
    $sysrcv set x [$targiv get x] 
    $sysrcv set y [$targiv get y] 
    $sysrcv set z [$targiv get z]

    # write the system-congruent output pdb
    [atomselect $sys all] writepdb $outpdb 
    
    set initrmsd [measure rmsd $a1 $a2]
    puts "Initial RMSD: [format "%0.2f" $initrmsd] Angstroms."
}

main $argv
exit
