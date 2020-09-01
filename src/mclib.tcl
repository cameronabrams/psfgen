# mclib
#
# library of custom routines for monte-carlo-based energy minimization 


# the actual pairwise interaction energy is calculated using a C-code
# swigged into a TcL procedure
load ${PSFGEN_BASEDIR}/src/energy.so energy

#Tcl functions for converting a Tcl double list to a C array
proc ListToArray {l} {
   set length [llength $l]
   set a [new_array $length]
   set i 0
   foreach item $l {
        array_setitem $a $i $item
        incr i 1
    }
  return $a
}

proc ListToArray_Data { a l } {
    set i 0
    foreach item $l {
	array_setitem $a $i $item
	incr i 1
    }
    return $a
}

#Tcl functions for converting a Tcl interger list to a C array
proc intListToArray {l} {
   set length [llength $l]
   set a [new_arrayint $length]
   set i 0
   foreach item $l {
        arrayint_setitem $a $i $item
        incr i 1
    }
  return $a
}

proc random {min max} {
    set sp [expr {$max - $min}]
    return [expr {rand() * $sp + $min}]
}

proc randomint {min max} {
    return [expr (int([random $min $max]))]
}

# read_charmm_nb
#
# reads all nonbonded pairwise interaction parameters (LJ)
# from a list of charmm-style parameter or toppar files
proc read_charmm_nb { par typ eps sig } {
  upvar $typ t
  upvar $eps e
  upvar $sig s

  foreach pfile $par {
    puts "Reading charmm file $pfile..."
    if { ! [file exists $pfile] } {
      puts "ERROR: $pfile does not exist"
      exit
    }
    set fp [open $pfile "r"]
    set data [read -nonewline $fp]
	
    close $fp
    set lines [split $data \n]

    set active 0
    foreach l $lines {
      set fc [string index $l 0]
      if { $fc != "!" } {
        if { ! [string compare [lindex $l 0] "HBOND"] } {
          set active 0
        }
        if { $active && [llength $l] } {
          set fw [lindex $l 0]
	  set fc2 [string index $fw 0]
	  if  { [string compare -nocase $fw "cutnb"] && $fc2 != "!" } {
	    if { ! [info exists t] } {
	      lappend t [lindex $l 0]
	      lappend e [lindex $l 2]
	      lappend s [lindex $l 3]
	      puts "adding typ [lindex $l 0] from $pfile..."
	    } elseif { [lsearch $t [lindex $l 0]] == -1 } {
              lappend t [lindex $l 0]
              lappend e [lindex $l 2]
	      lappend s [lindex $l 3]
	      puts "adding typ [lindex $l 0] from $pfile..."
	    }
	  }
        }
        if { ![string compare [lindex $l 0] "NONBONDED"] } {
          set active 1
        }
      }
    } 
  }
}

# initialize_dataspace
#
# sets up an initial dataspace for transferring coordinates to/from TcL and C
proc initialize_dataspace { mmr lgr } {
    set dataspace {}
    lappend dataspace [ListToArray [$mmr get x]]
    lappend dataspace [ListToArray [$mmr get y]]
    lappend dataspace [ListToArray [$mmr get z]]
    lappend dataspace [ListToArray [$lgr get x]]
    lappend dataspace [ListToArray [$lgr get y]]
    lappend dataspace [ListToArray [$lgr get z]]
    return $dataspace
}

# UNDER CONSTRUCTION
# my_nonbonded_interaction_energy
#
# computes and returns the nonbonded interaction energy between atoms in the "mmr" and "lgr" 
# atomselections.  "roughenergy" is a C function defined in "energy.c".
proc my_nonbonded_interaction_energy { mmr lgr epsilon sigma index1 index2 ch1 ch2 sizemmr sizelgr dataspace } {
    
    set E 0.0

    ListToArray_Data [lindex $dataspace 0] [$mmr get x]
    ListToArray_Data [lindex $dataspace 1] [$mmr get y]
    ListToArray_Data [lindex $dataspace 2] [$mmr get z]
    ListToArray_Data [lindex $dataspace 3] [$lgr get x]
    ListToArray_Data [lindex $dataspace 4] [$lgr get y]
    ListToArray_Data [lindex $dataspace 5] [$lgr get z]

    set E [roughenergy [lindex $dataspace 0] [lindex $dataspace 1] [lindex $dataspace 2] [lindex $dataspace 3] [lindex $dataspace 4] [lindex $dataspace 5] $epsilon $sigma $index1 $index2 $ch1 $ch2 $sizemmr $sizelgr]

    return $E
}


proc restr_make_selections { restr lig_id sub_id } {
    foreach r $restr {
	set ls [atomselect $lig_id "[lindex $r 0]"]
	$ls global
	set ss [atomselect $sub_id "[lindex $r 1]"]
	$ss global
	lappend rr [list $ls $ss [lindex $r 2] [lindex $r 3]]
    }
    return $rr
}

proc measure_restr { restr m prefix } {
    set s1 {}
    set s2 {}
    set pairs {}
    set fp [open "${prefix}.dat" "w"]
    puts -nonewline $fp "\# "
    set n 0
    foreach r $restr {
#	puts "DB: restr $n : $r : [lindex $r 0]"
	set sel1 [atomselect $m "[lindex $r 0]"]
	lappend s1 $sel1
	set sel2 [atomselect $m "[lindex $r 1]"]
	lappend s2 $sel2
#	puts "DB: restr $n ligand $sel1 [$sel1 num] substrate $sel2 [$sel2 num]"
	puts -nonewline $fp "[$sel1 get resname][$sel1 get resid][$sel1 get name]-[$sel2 get resname][$sel2 get resid][$sel2 get name] "
	set pair [list [$sel1 get index] [$sel2 get index]]
	lappend pairs $pair
	incr n
    }
#    puts "DB: pairs : $pairs"
    puts $fp ""
    set data {}
    foreach p $pairs {
	set td [measure bond $p frame all]
	puts "[llength $td] data."
	lappend data $td
    }

    for {set i 0} { $i < [llength [lindex $data 0]] } { incr i } {
	puts -nonewline $fp "$i "
	foreach d $data {
	    puts -nonewline $fp "[lindex $d $i] "
	}
	puts $fp ""
    }
    close $fp
}

proc measure_intermol_restr { restr m1 m2 frame } {
    set s1 {}
    set s2 {}
    set R {}
    foreach r $restr {
	set sel1 [atomselect $m1 "[lindex $r 0]" frame $frame]
	lappend s1 $sel1
	set sel2 [atomselect $m2 "[lindex $r 1]" frame $frame]
	lappend s2 $sel2
	set dr {}
	lappend dr [expr [$sel1 get x] - [$sel2 get x]]
	lappend dr [expr [$sel1 get y] - [$sel2 get y]]
	lappend dr [expr [$sel1 get z] - [$sel2 get z]]
	lappend R [veclength $dr]
    }
    return $R
}

proc get_wt { val min max } {
    return [expr ($val-$min)/($max-$min)]
}


proc read_trajectory { traj } {
    set psf [lindex $traj 0]
    set coord [lindex $traj 1]
    set beg [lindex $traj 2]
    mol new $psf
    mol addfile $coord waitfor all autobonds off first $beg
    return [molinfo top get id]
}

proc restraint_energy { cyc restr sub_frame lig_frame} {
    set E 0.0
    foreach r $restr {
	set la [lindex $r 0]
	set sa [lindex $r 1]

	$la frame $lig_frame
	$sa frame $sub_frame

	set sa_r [list [$sa get x] [$sa get y] [$sa get z]]
	set la_r [list [$la get x] [$la get y] [$la get z]]

	set R [veclength [vecsub $sa_r $la_r]]
	set eq [lindex $r 2]
	set k [lindex $r 3]
	set thisE [expr  $k*pow($R-$eq,2)]
	set E [expr $E + $thisE]

    }
    return $E
}

proc get_move { x move_weights } {
    set i 0
    set j -1
    foreach m $move_weights {
	set wt [lindex $m 1]
	if {$j == -1 && [expr $x < $wt]} {
	    set j $i
	}
	incr i
    }
    if {$j > -1} {
	return [lindex [lindex $move_weights $j] 0]
    } else {
	return INVALID
    }
}

proc mc_read_init { fn Seed Par Ligtraj Subtraj Ligdocksel Subdocksel Subplacesel Ligfirstframe Subfirstframe Restr Cycles KT DAngleDeg Rmin Rmax Move_weights OutputInterval } {
    upvar $Seed seed
    upvar $Par par
    upvar $Ligtraj ligtraj
    upvar $Subtraj subtraj
    upvar $Ligdocksel lds
    upvar $Subdocksel sds
    upvar $Subplacesel sps
    upvar $Ligfirstframe ligfirstframe
    upvar $Subfirstframe subfirstframe
    upvar $Restr restr
    upvar $Cycles cycles
    upvar $KT kT
    upvar $DAngleDeg dAngleDeg
    upvar $Rmin rmin
    upvar $Rmax rmax
    upvar $Move_weights move_weights
    upvar $OutputInterval oi

    set restr {}

    set fp [open "$fn" "r"]
    set data [read -nonewline $fp]
    set lines [split $data \n]
    foreach l $lines {
#	puts "DB: elements [llength $l] in line: $l"
	set key [lindex $l 0]
	set cmt [string first \# $key]
	if {$cmt == 0} {
	    # do nothing; it's a comment
	} else {
	    switch $key {
		SEED {
		    set seed [lindex $l 1]
		    if { $seed == "*" } {
			set seed [clock clicks]
		    }
		}
		PAR {
		    lappend par [lindex $l 1]
		}
		LIGTRAJ {
		    set ligtraj [list [lindex $l 1] [lindex $l 2] [lindex $l 3]]
		}
		LIGFIRSTFRAME {
		    set ligfirstframe [lindex $l 1]
		}
		SUBTRAJ {
		    set subtraj [list [lindex $l 1] [lindex $l 2] [lindex $l 3]]
		}
		SUBFIRSTFRAME {
		    set subfirstframe [lindex $l 1]
		}
		LIGDOCKSEL {
		    set lds [lindex $l 1]
		}
		SUBDOCKSEL {
		    set sds [lindex $l 1]
		}
		SUBPLACESEL {
		    set sps [lindex $l 1]
		}
		RESTR {
		    lappend restr [list [lindex $l 1] [lindex $l 2] [lindex $l 3] [lindex $l 4]]
		}
		CYCLES {
		    set cycles [lreplace $l 0 0]
		}
		KT {
		    set kT [lindex $l 1]
		}
		DANGLEDEG {
		    set dAngleDeg [lindex $l 1]
		}
		RMIN {
		    set rmin [lindex $l 1]
		}
		RMAX {
		    set rmax [lindex $l 1]
		}
		CWEIGHT {
		    lappend move_weights [list [lindex $l 1] [lindex $l 2]]
		}
		OUTPUTINTERVAL {
		    set oi [lindex $l 1]
		}
		default {
		    puts "not ready to handle key $key"
		}
	    }
	}
    }
}

proc mc_report_init { seed par lt st lds sds sps ligfirstframe subfirstframe restr cycles kT dAngleDeg rmin rmax move_weights outputInterval } {
    puts "INFO: seed: $seed"
    puts "INFO: par: $par"
    puts "INFO: ligtraj: $lt"
    puts "INFO: ligfirstframe: $ligfirstframe"
    puts "INFO: subtraj: $st"
    puts "INOF: subfirstframe: $subfirstframe"
    puts "INFO: ligand docking selection: $lds"
    puts "INFO: substrate docking selection: $sds"
    puts "INFO: substrate initial placement selection: $sps"
    if {[llength $restr]} {
	set n 0
	foreach r $restr {
	    puts "INFO: restr $n: ligand_atom ([lindex $r 0]) substrate_atom ([lindex $r 1]) L0 [lindex $r 2] k [lindex $r 3]"
	    incr n
	}
    }
    puts "INFO: cycles: $cycles"
    puts "INFO: kT: $kT"
    puts "INFO: dAngleDeg: $dAngleDeg"
    puts "INFO: rmin/rmax: ${rmin}/${rmax}"
    set n 0
    foreach m $move_weights {
	puts "INFO: move $n weight [lindex $m 0] [lindex $m 1]"
	incr n
    }
    puts "INFO: output interval: $outputInterval cycles"
}

proc initialize_ligand_position { lig_id sub_id lig_sel sub_sel } {
    global M_PI
    
    set s1 [atomselect $sub_id "name CA and ($sub_sel)"]
    set cm [measure center $s1]
    set s2 [atomselect $lig_id "$lig_sel"]
    $s1 frame 0
    $s2 frame 0
    set lig_nf [molinfo $lig_id get numframes]
    set z1 [random -180.0 180.0]
    set z2 [random -180.0 180.0]
    set cx [random -1 1]
    set x [expr acos($cx)]
    set z [random 0 1]
    if { [expr $z < 0.5 ] } {
	set x [expr -$x]
    }
    set x [expr {$x * 57.296}]
    set Eulerz1 [trans center $cm axis z $z1 deg]
    set Eulerx  [trans center $cm axis x $x deg]
    set Eulerz2 [trans center $cm axis z $z2 deg]
    set R [transmult $Eulerz2 $Eulerx $Eulerz1]

    puts "INFO: Initializing ligand position:"
    puts "INFO: Ligand (molid $lig_id frames $lig_nf)  : $lig_sel "
    puts "INFO: Raw dock selection on substrate ($sub_id) : $sub_sel"
    puts "INFO: Raw dock selection center at $cm"
    for {set i 0 } { $i < $lig_nf} {incr i} {
	$s2 frame $i
	$s2 moveby [vecsub $cm [measure center $s2]]
	$s2 move $R
    }    
}

