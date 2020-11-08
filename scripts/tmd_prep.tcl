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
# In addition to identifying the three input files (PSF, PDB or COOR, and PDB),
# the user must specify the alignment basis and the target selection.  
#
# ALIGNMENT BASIS.  This is a comma-separated list of colon-separated pairs.  
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
# Then their "occupancy" flags are set to 1.
#
# A full selection of all atoms in the system is then output as the
# TMD PDB.
#
# Cameron Abrams cfa22@drexel.edu

# generates an atomselect string from the command-line shortcode
# shortcode format is ABC_xxx-yyy;zzz-ggg;....
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
    # set stringmap [array get arr]
    # puts "stringmap $stringmap"
#     puts "shortcode $shortcode"
    # set newshortcode [string map $stringmap $shortcode]
    # puts "newshortcode $newshortcode"
     set myl [split $shortcode " "]
#     puts "myl"
     for {set i 0} {$i<[llength $myl]} {incr i} {
	set rawwrd [lindex $myl $i]
	set newwrd $rawwrd
	set chrs [split $rawwrd ""]
        set anum ""
	set fst -1
	set lst 0
	for {set j 0} {$j<[llength $chrs]} {incr j} {
	    set c [lindex $chrs $j]
#	    puts -nonewline "-$c"
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
#	    puts -nonewline "(anum $anum map $arr($anum))"
        }
#        puts -nonewline "|"
     }
#     puts ""
     set newshortcode [join $myl " "]
#     puts "newshortcode $newshortcode"
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
    mol new $tarpdb 
    set tar [molinfo top get id]
    set alb [build_atom_select_strings $albs]
    set tas [build_atom_select_strings $tass]
    # if the system and target have different resid numbering, it is necessary to 
    # map the target resids to their counterpart values in the system; we
    # refer to this as a "21" mapping; "2" indicates target.  This involves
    # changing resid values in the target pdb molecule and changing them
    # in the atomselections against that molecule.
    # For example, if the system is based on 5x58 and the target is 6vyb
    # we will renumber the residues in the 6vyb molecule instance by asking
    # what resid in 5x58 corresponds to a particular resid in 6vyb, or a 
    # "6vyb-to-5x58" map.
    if { [info exists mymap ] } {
	# renumber resid's in the target pdb molecule
	set newresid [list]
	set tarall [atomselect $tar "protein"]
	set oldresid [$tarall get resid]
	foreach r $oldresid {
	   lappend newresid $mymap($r)
	}
	puts "$newresid"
	$tarall set resid $newresid
	# change atomselection resids in the target alignment basis
	set mapped [mapresids [lindex $alb 1] mymap]
	lset alb 1 $mapped
	# change atomselection resids in the target selection
	set mapped [mapresids [lindex $tas 1] mymap]
	lset tas 1 $mapped
    }
    set tarall [atomselect $tar all]
    set taralb [atomselect $tar "$seland and ([lindex $alb 1])"]
    set sysalb [atomselect $sys "$seland and ([lindex $alb 0])"]
    puts "Alignment basis: SYS([lindex $alb 0]) and TARG([lindex $alb 1])"
    puts "Target selection: SYS([lindex $tas 0]) and TARG([lindex $tas 1])"
    puts "System: $psf $cor"
    puts "Target pdb: $tarpdb"
    puts "taralb num [$taralb num] sysalb num [$sysalb num]"
    $tarall move [measure fit $taralb $sysalb]
    set sysrcv [atomselect $sys "$seland and ([lindex $tas 0])"]
    set targiv [atomselect $tar "$seland and ([lindex $tas 1])"]
    puts "targiv num [$targiv num] sysrcv num [$sysrcv num]"
    $sysrcv set x [$targiv get x] 
    $sysrcv set y [$targiv get y] 
    $sysrcv set z [$targiv get z]
    $sysrcv set occupancy 1
    [atomselect $sys all] writepdb $outpdb 
     
}

main $argv
exit
