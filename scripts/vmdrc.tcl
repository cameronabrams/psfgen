# .vmdrc
# cameron f abrams, 2018
# cfa22@drexel.edu

# draw an arrow from pt 'start' to point 'end'
proc vmd_draw_arrow {mol start end} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
    graphics $mol cylinder $start $middle radius 0.15
    graphics $mol cone $middle $end radius 0.25
}

# draw a sphere at the selection's center of mass
proc showcom { sel mol } {
    set nf [molinfo $mol get numframes]
    graphics $mol color red
    for {set i 0} {$i < $nf} {incr i} {
	$sel frame $i
	graphics $mol sphere [measure center $sel weight mass]
    }
}

# move a selection by vector 'offset'
proc moveby {sel offset} {
  foreach coord [$sel get {x y z}] {
    lappend newcoords [vecadd $coord $offset]
  }
  $sel set {x y z} $newcoords
}

# calculate and report box dimensions based on top molecule
proc getpbc { pad } {
    if {![info exists pad]} {
	set pad 0.1
    }
    set sel [atomselect top all]
    $sel frame 0
    set mm [measure minmax $sel]
    set ll [lindex $mm 0]
    set ur [lindex $mm 1]
    set sp [vecsub $ur $ll]
    puts "cellBasisVector1   [expr [lindex $sp 0] + $pad] 0.0 0.0"
    puts "cellBasisVector2   0.0 [expr [lindex $sp 1] + $pad] 0.0"
    puts "cellBasisVector3   0.0 0.0 [expr [lindex $sp 2] + $pad]"
    puts "cellOrigin         [measure center $sel]"
    return "cellBasisVector1   [expr [lindex $sp 0] + $pad] 0.0 0.0\ncellBasisVector2   0.0 [expr [lindex $sp 1] + $pad] 0.0\ncellBasisVector3   0.0 0.0 [expr [lindex $sp 2] + $pad]\ncellOrigin         [measure center $sel]"
}

# compute and return internal coordinate (IC) for the four atoms listed (by index)
proc getic { at1 at2 at3 at4 } {
    set b12 [measure bond [list $at1 $at2]]
    set b23 [measure bond [list $at2 $at3]]
    set b34 [measure bond [list $at3 $at4]]
    set a123 [measure angle [list $at1 $at2 $at3]]
    set a234 [measure angle [list $at2 $at3 $at4]]
    set d1234 [measure dihed [list $at1 $at2 $at3 $at4]]
    puts "[format %.3f $b12] [format %.3f $a123] [format %.3f $d1234] [format %.3f $a234] [format %.3f $b34]"
}

proc geticimp { at1 at2 at3 at4 } {
    set b13 [measure bond [list $at1 $at3]]
    set b23 [measure bond [list $at2 $at3]]
    set b34 [measure bond [list $at3 $at4]]
    set a132 [measure angle [list $at1 $at3 $at2]]
    set a234 [measure angle [list $at2 $at3 $at4]]
    set d1234 [measure dihed [list $at1 $at2 $at3 $at4]]
    puts "[format %.3f $b13] [format %.3f $a132] [format %.3f $d1234] [format %.3f $a234] [format %.3f $b34]"
}

# compute and return IC for four atoms listed by name
proc geticN { an1 an2 an3 an4 } {
    set as1 [atomselect top "name $an1"]
    set at1 [$as1 get index]
    $as1 delete
    set as2 [atomselect top "name $an2"]
    set at2 [$as2 get index]
    $as2 delete
    set as3 [atomselect top "name $an3"]
    set at3 [$as3 get index]
    $as3 delete
    set as4 [atomselect top "name $an4"]
    set at4 [$as4 get index]
    $as4 delete
    puts -nonewline "IC $an1 $an2 $an3 $an4 "
    getic $at1 $at2 $at3 $at4
}

proc geticimpN { an1 an2 an3 an4 } {
    set as1 [atomselect top "name $an1"]
    set at1 [$as1 get index]
    $as1 delete
    set as2 [atomselect top "name $an2"]
    set at2 [$as2 get index]
    $as2 delete
    set as3 [atomselect top "name $an3"]
    set at3 [$as3 get index]
    $as3 delete
    set as4 [atomselect top "name $an4"]
    set at4 [$as4 get index]
    $as4 delete
    puts -nonewline "IC $an1 $an2 *$an3 $an4 "
    geticimp $at1 $at2 $at3 $at4
}

# detect pucker of all proline residues in list Resids on chain 'chain' of molecule 'molid'
proc propuck { molid chain Resids } {
    set result {}
    foreach resid $Resids {
	set c  [atomselect $molid "chain $chain and resid $resid and name C"]
	set ca [atomselect $molid "chain $chain and resid $resid and name CA"]
	set cb [atomselect $molid "chain $chain and resid $resid and name CB"]
	set cg [atomselect $molid "chain $chain and resid $resid and name CG"]
	set cd [atomselect $molid "chain $chain and resid $resid and name CD"]
	set n  [atomselect $molid "chain $chain and resid $resid and name N"]
	set hb1 [atomselect $molid "chain $chain and resid $resid and name HB1"]
	set cn [atomselect $molid "chain $chain and resid [expr $resid + 1] and name C"]

	set x1 [measure dihed [list [$c get index] [$ca get index] [$cb get index] [$cg get index]]]
	set x2 [measure dihed [list [$ca get index] [$cb get index] [$cg get index] [$cd get index]]]
	set x3 [measure dihed [list [$cb get index] [$cg get index] [$cd get index] [$n get index]]]
	set x4 [measure dihed [list [$cg get index] [$cd get index] [$n get index] [$cn get index]]]

	set t1 [measure angle [list [$hb1 get index] [$cb get index] [$cg get index]]]

	set res DOWN
	# criterion of Milner-White, Bell, and Maccallum, J. Mol. Biol. 228:725-734 (1992)
	set critv [expr $x1 + $x3 - $x2 -$x4]
	if {[expr $critv > -100.0]} {
	    set res UP
	}
#	puts "DB: [$ca get resname][$ca get resid]: x1 [format %.2f $x1] x2 [format %.2f $x2] x3 [format %.2f $x3] x4 [format %.2f $x4] critv [format %.2f $critv] : $res : t1 [format %2f $t1]"
	lappend result $res
    }
    return $result
}

# draw a cage from corner ll to corner ur
proc draw_cage_cc {mol radius ll ur color} {
    set p000 [list [lindex $ll 0] [lindex $ll 1] [lindex $ll 2]]
    set p100 [list [lindex $ur 0] [lindex $ll 1] [lindex $ll 2]]
    set p110 [list [lindex $ur 0] [lindex $ur 1] [lindex $ll 2]]
    set p010 [list [lindex $ll 0] [lindex $ur 1] [lindex $ll 2]]
    set p001 [list [lindex $ll 0] [lindex $ll 1] [lindex $ur 2]]
    set p101 [list [lindex $ur 0] [lindex $ll 1] [lindex $ur 2]]
    set p111 [list [lindex $ur 0] [lindex $ur 1] [lindex $ur 2]]
    set p011 [list [lindex $ll 0] [lindex $ur 1] [lindex $ur 2]]

    graphics $mol color $color
    graphics $mol sphere $p000 radius $radius
    graphics $mol sphere $p100 radius $radius
    graphics $mol sphere $p110 radius $radius
    graphics $mol sphere $p010 radius $radius
    graphics $mol sphere $p001 radius $radius
    graphics $mol sphere $p101 radius $radius
    graphics $mol sphere $p111 radius $radius
    graphics $mol sphere $p011 radius $radius

    graphics $mol cylinder $p000 $p100 radius $radius
    graphics $mol cylinder $p100 $p110 radius $radius
    graphics $mol cylinder $p110 $p010 radius $radius
    graphics $mol cylinder $p010 $p000 radius $radius
    graphics $mol cylinder $p001 $p101 radius $radius
    graphics $mol cylinder $p101 $p111 radius $radius
    graphics $mol cylinder $p111 $p011 radius $radius
    graphics $mol cylinder $p011 $p001 radius $radius
    graphics $mol cylinder $p000 $p001 radius $radius
    graphics $mol cylinder $p100 $p101 radius $radius
    graphics $mol cylinder $p110 $p111 radius $radius
    graphics $mol cylinder $p010 $p011 radius $radius

    
}

# draw a cage around a molecule
proc draw_cage {mol radius} {

    graphics $mol delete all

    set all [atomselect $mol "all"]

    if {1} {
	set cc [molinfo $mol get center]
	set center [lindex $cc 0]
	set ao [lindex $center 0]
	set bo [lindex $center 1]
	set co [lindex $center 2]
	set a [molinfo $mol get a]
	set b [molinfo $mol get b]
	set c [molinfo $mol get c]
	set a2 [expr $a/2]
	set b2 [expr $b/2]
	set c2 [expr $c/2]
	set ll [list [expr $ao-$a2] [expr $bo-$b2] [expr $co-$c2]]
	set ur [list [expr $ao+$a2] [expr $bo+$b2] [expr $co+$c2]]
    } else {
	set minmax [measure minmax $all]
	set ll [lindex $minmax 0]
	set ur [lindex $minmax 1]
    }

    set sp [vecsub $ur $ll]
    puts [format "%.2f %.2f %.2f" [lindex $sp 0] [lindex $sp 1] [lindex $sp 2]]

    set p000 [list [lindex $ll 0] [lindex $ll 1] [lindex $ll 2]]
    set p100 [list [lindex $ur 0] [lindex $ll 1] [lindex $ll 2]]
    set p110 [list [lindex $ur 0] [lindex $ur 1] [lindex $ll 2]]
    set p010 [list [lindex $ll 0] [lindex $ur 1] [lindex $ll 2]]
    set p001 [list [lindex $ll 0] [lindex $ll 1] [lindex $ur 2]]
    set p101 [list [lindex $ur 0] [lindex $ll 1] [lindex $ur 2]]
    set p111 [list [lindex $ur 0] [lindex $ur 1] [lindex $ur 2]]
    set p011 [list [lindex $ll 0] [lindex $ur 1] [lindex $ur 2]]


    graphics $mol color 4
    graphics $mol sphere $p000 radius $radius
    graphics $mol sphere $p100 radius $radius
    graphics $mol sphere $p110 radius $radius
    graphics $mol sphere $p010 radius $radius
    graphics $mol sphere $p001 radius $radius
    graphics $mol sphere $p101 radius $radius
    graphics $mol sphere $p111 radius $radius
    graphics $mol sphere $p011 radius $radius

    graphics $mol cylinder $p000 $p100 radius $radius
    graphics $mol cylinder $p100 $p110 radius $radius
    graphics $mol cylinder $p110 $p010 radius $radius
    graphics $mol cylinder $p010 $p000 radius $radius
    graphics $mol cylinder $p001 $p101 radius $radius
    graphics $mol cylinder $p101 $p111 radius $radius
    graphics $mol cylinder $p111 $p011 radius $radius
    graphics $mol cylinder $p011 $p001 radius $radius
    graphics $mol cylinder $p000 $p001 radius $radius
    graphics $mol cylinder $p100 $p101 radius $radius
    graphics $mol cylinder $p110 $p111 radius $radius
    graphics $mol cylinder $p010 $p011 radius $radius


}

# shift center of mass
proc acm { molid } {
    set nf [molinfo top get numframes]
    set p [atomselect $molid protein]
    for { set i 0 } { $i < $nf } { incr i } {
	$p frame $i
	$p moveby [vecscale -1 [measure center $p weight mass]]
    }
    $p delete
}

# align whole molecule
proc awm { molid } {
    set nf [molinfo top get numframes]
    set p [atomselect $molid "protein or (ion and within 4.0 of protein)"]
    set pr [atomselect $molid "protein or (ion and within 4.0 of protein)"]
    $pr frame 0
    for { set i 0 } { $i < $nf } { incr i } {
	$p frame $i
	$p move [measure fit $p $pr]
    }
    $p delete
    $pr delete
}


# align molecule on selection
proc amos { molid selstr } {
    set nf [molinfo top get numframes]
    set p [atomselect $molid protein]
    set sr [atomselect $molid "$selstr"]
    set s [atomselect $molid "$selstr"]
    $sr frame 0
    for { set i 0 } { $i < $nf } { incr i } {
	$p frame $i
	$s frame $i
	$p move [measure fit $s $sr]
    }
    $p delete
    $s delete
    $sr delete
}

# align molecule on selection
proc amoss { molid selstr sel2str } {
    set nf [molinfo top get numframes]
    set p [atomselect $molid "$sel2str"]
    set sr [atomselect $molid "$selstr"]
    set s [atomselect $molid "$selstr"]
    $sr frame 0
    for { set i 0 } { $i < $nf } { incr i } {
	$p frame $i
	$s frame $i
	$p move [measure fit $s $sr]
    }
    $p delete
    $s delete
    $sr delete
}

proc fwtaln { molid } {
    set nf [molinfo top get numframes]
    set p [atomselect $molid protein]
    set q [atomselect $molid protein]
    for { set i 1 } { $i < $nf } { incr i } {
	set rf [expr $i - 1]
	$p frame $rf
	$q frame $i
	set tm [measure fit $q $p]
	$q move $tm
    }
    $p delete
    $q delete
}

# this is really dumb:
proc bwtaln { molid } {
    set nf [molinfo top get numframes]
    set p [atomselect $molid protein]
    set q [atomselect $molid protein]
    for { set i [expr $nf - 2] } { $i >= 0 } { incr i -1 } {
	puts "$i / $nf"
	set af [expr $i + 1]
	$p frame $i
	$q frame $af
	set tm [measure fit $q $p]
	for { set j $af } { $j < $nf } { incr j } {
	    puts " ... $j / $nf"
	    $q frame $j
	    $q move $tm
	}
    }
    $p delete
    $q delete
}

proc get_azim_angle { x y } {
    set angle [expr atan2($y,$x)]

    if {[expr $y < 0]} {
	set angle [expr $angle + 2 * 3.1415928]
    }
    return $angle
}

# this routine aligns principle axes by direct rotations
proc alnpa { molid } {
    set nf [molinfo $molid get numframes]
    set a [atomselect $molid protein]
    for {set f 0} {$f < $nf} {incr f} {
	$a frame $f
	set dd [measure inertia $a]
	puts "$dd"
	set c [lindex $dd 0]
	set m [lindex $dd 1]
	set x0 [lindex $m 2] ; # longest axis
	puts "$c $x0"
	# move to center
	$a moveby [vecscale $c -1]
	# rotate on global z to bring intertial z into global z-x plane
	set x [lindex $x0 0]
	set y [lindex $x0 1]
	set zangle [expr -180.0/3.1415928 * [get_azim_angle $x $y]]
	puts "z1 rot by $zangle deg..."
	$a move [trans center {0 0 0} axis z $zangle deg]
	
	# convention, let longest axis be z, then adopt a rh coordinate system to pick x and y
	set dd [measure inertia $a]
	set c [lindex $dd 0]
	set m [lindex $dd 1]
	set x0 [lindex $m 2] ; # longest axis
	
	set yangle [expr -180.0/3.1415928 * [expr acos([lindex $x0 2])]]
	puts "y rot by $yangle deg..."
	$a move [trans center {0 0 0} axis y $yangle deg]

	set dd [measure inertia $a]
	set c [lindex $dd 0]
	set m [lindex $dd 1]
	set x0 [lindex $m 2] ; # longest axis
	set x1 [lindex $m 1] ; # 
	set x2 [lindex $m 0] ; # shortest axis
	set c1 [veccross $x1 $x2]
	set d1 [vecdot $c1 $x0]
	if {[expr $d1 > 0]} {
	    set local_x $x1
	    set local_y $x2
	} else {
	    set local_x $x2
	    set local_y $x1
	}
	
	set x [lindex $local_x 0]
	set y [lindex $local_x 1]
	set zangle [expr -180.0/3.1415928 * [get_azim_angle $x $y]]
	puts "z2 rot by $zangle deg.."
	$a move [trans center {0 0 0} axis z [expr -1 * $zangle] deg]
	puts "$f done."
	
    }
}

# secondary structure measurement/rendering
# start the cache for a given molecule
proc start_sscache { molid } {
    global sscache_data
    global vmd_frame
    # set a trace to detect when an animation frame changes
    trace variable vmd_frame($molid) w sscache_from_trace
    return
}

proc stop_sscache { molid } { 
    global vmd_frame
    puts "stopping sscache at frame [molinfo $molid get frame] molid $molid..."
    trace vdelete vmd_frame($molid) w sscache_from_trace
    return
}

# reset the whole secondary structure data cache
proc reset_sscache {} {
    if [info exists sscache_data] {
	unset sscache_data
    }
    return
}

# when the frame changes, trace calls this function
proc sscache_from_trace { name1 name2 op } {
    global sscache_data

    # undo the tcl trace argument naming convention
    set molid $name2
    set frame [molinfo $molid get frame]

    # get the protein CA atoms
    set sel [atomselect $molid "protein name CA"]
    
    # see if the ss data exists in the cache
    if [info exists sscache_data($molid,$frame)] {
	$sel set structure $sscache_data($molid,$frame)
	return
    }
	
    # doesn't exist, so (re)calculate it
    vmd_calculate_structure $molid
    # save the data for next time
    set sscache_data($molid,$frame) [$sel get structure]

    return
}	

proc sstrace_all { molid filename } {
    global sscache_data

    set nf [molinfo $molid get numframes]
    set fp [open $filename "w"]
    
    for {set f 0} {$f < $nf} {incr f} {
        foreach L $sscache_data($molid,$f) {
            #$puts -nonewline $fp "$f "
            switch $L {
                H {
                    set v 4
                }
                E {
                    set v 3
                }
                T {
                    set v 2
                }
                B {
                    set v 1
                }
                C {
                    set v 0
                }
                default {
                    set v 0
                }
            }
            puts -nonewline $fp "$v "
        }
        puts $fp ""
    }
    close $fp
}

proc sstrace_line { molid filename } {
    global sscache_data

    set nf [molinfo $molid get numframes]
    set fp [open $filename "w"]
    set sel [atomselect $molid "protein name CA"]
    for {set f 0} {$f < $nf} {incr f} {
        set r 0
        # see if the ss data exists in the cache
        if [info exists sscache_data($molid,$f)] {
            $sel set structure $sscache_data($molid,$f)
        } else {
            # doesn't exist, so (re)calculate it
            animate goto $f
            vmd_calculate_structure $molid
            # save the data for next time
            set sscache_data($molid,$f) [$sel get structure]
        }
        puts -nonewline $fp "$f "
        foreach L $sscache_data($molid,$f) {
            puts -nonewline $fp "$L"
        }
        puts $fp ""
    }
    close $fp
    puts "Created $filename."
}

proc sstrace_gnuplot { molid filename force_recalc } {
    global sscache_data

    set nf [molinfo $molid get numframes]
    set fp [open $filename "w"]
    set sel [atomselect $molid "protein name CA"]
    for {set f 0} {$f < $nf} {incr f} {
        set r 0
        # see if the ss data exists in the cache
        if {!$force_recalc && [info exists sscache_data($molid,$f)]} {
            $sel set structure $sscache_data($molid,$f)
        } else {
            # doesn't exist, so (re)calculate it
            animate goto $f
            vmd_calculate_structure $molid
            # save the data for next time
            set sscache_data($molid,$f) [$sel get structure]
        }
        foreach L $sscache_data($molid,$f) {
            puts -nonewline $fp "$f $r "
            switch $L {
                H {
                    set v 2
                }
                E {
                    set v -2
                }
                T {
                    set v -1
                }
                B {
                    set v 1
                }
                C {
                    set v 0
                }
                default {
                    set v 0
                }
            }
            puts $fp "$v"
            incr r
        }
        puts $fp ""
    }
    close $fp
    puts "Created $filename."
}

# move whole trajectories
proc tmoveby { molid sel vec } {
    set nf [molinfo $molid get numframes]
    for {set i 0} {$i < $nf} {incr i} {
	$sel frame $i
	$sel moveby $vec 
    }
}

proc tmove { molid sel mat } {
    set nf [molinfo $molid get numframes]
    for {set i 0} {$i < $nf} {incr i} {
	$sel frame $i
	$sel move $mat
    }
}

# compute vector that will rotate current view around screen z-axis
proc vzv { molid } {
    set eye_vector [vectrans \
			[measure inverse [lindex [molinfo $molid get rotate_matrix] 0]] \
			{0 0 1}
		   ]
    return $eye_vector
}

# compute vector that will rotate current view around screen x-axis
proc vxv { molid } {
    set eye_vector [vectrans \
			[measure inverse [lindex [molinfo $molid get rotate_matrix] 0]] \
			{1 0 0}
		   ]
    return $eye_vector
}

# compute vector that will rotate current view around screen y-axis
proc vyv { molid } {
    set eye_vector [vectrans \
			[measure inverse [lindex [molinfo $molid get rotate_matrix] 0]] \
			{0 1 0}
		   ]
    return $eye_vector
}

# return 1-letter amino acid code given the 3-character code
proc aa_321 { aa3 } {
    set aa1 {}
    foreach aa $aa3 {
        switch $aa {
            ALA {
                set a1 A
            }
            VAL {
                set a1 V
            }
            ILE {
                set a1 I
            }
            LEU {
                set a1 L
            }
            PRO {
                set a1 P
            }
            PHE {
                set a1 F
            }
            TRP {
                set a1 W
            }
            MET {
                set a1 M
            }
            GLY {
                set a1 G
            }
            SER {
                set a1 S
            }
            THR {
                set a1 T
            }
            CYS {
                set a1 C
            }
            ASN {
                set a1 N
            }
            GLN {
                set a1 Q
            }
            TYR {
                set a1 Y
            }
            HIS -
            HSD -
            HSE -
            HSP {
                set a1 H
            }
            LYS {
                set a1 K
            }
            ARG {
                set a1 R
            }
            ASP {
                set a1 D
            }
            GLU {
                set a1 E
            }
            default {
		puts "aa? $aa"
                set a1 "?"
            }
        }
        lappend aa1 $a1
    }
    return $aa1
}

atomselect macro dppc_head "resname DPPC and name C1 HA HB C11 H11A H11B C12 H12A H12B C13 H13A H13B H13C C14 H14A H14B H14C C15 H15A H15B H15C P O11 O12 O13 O14 N"
atomselect macro dppc_tail "resname DPPC and not name C1 HA HB C11 H11A H11B C12 H12A H12B C13 H13A H13B H13C C14 H14A H14B H14C C15 H15A H15B H15C P O11 O12 O13 O14 N"
atomselect macro glycan "resname NAG MAN BMA FUC GAL BGNA AMAN BMAN AFUC BGAL" 
menu main on

