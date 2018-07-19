proc vmd_draw_arrow {mol start end} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
    graphics $mol cylinder $start $middle radius 0.15
    graphics $mol cone $middle $end radius 0.25
}

proc ied {} {
        gopython /usr/local/lib/vmd/plugins/noarch/python/ied/ied.py
}

proc center_of_mass {selection} {
        # some error checking
        if {[$selection num] <= 0} {
                error "center_of_mass: needs a selection with atoms"
        }
        # set the center of mass to 0
        set com [veczero]
        # set the total mass to 0
        set mass 0
        # [$selection get {x y z}] returns the coordinates {x y z} 
        # [$selection get {mass}] returns the masses
        # so the following says "for each pair of {coordinates} and masses,
        #  do the computation ..."
        foreach coord [$selection get {x y z}] m [$selection get mass] {
           # sum of the masses
           set mass [expr $mass + $m]
           # sum up the product of mass and coordinate
           set com [vecadd $com [vecscale $m $coord]]
        }
        # and scale by the inverse of the number of atoms
        if {$mass == 0} {
                error "center_of_mass: total mass is zero"
        }
        # The "1.0" can't be "1", since otherwise integer division is done
        return [vecscale [expr 1.0/$mass] $com]
}

proc showcom { sel mol } {
    set nf [molinfo $mol get numframes]
    graphics $mol color red
    for {set i 0} {$i < $nf} {incr i} {
	$sel frame $i
	graphics $mol sphere [measure center $sel weight mass]
    }
}

proc moveby {sel offset} {
  foreach coord [$sel get {x y z}] {
    lappend newcoords [vecadd $coord $offset]
  }
  $sel set {x y z} $newcoords
}

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


if {[file exists /home/cfa/tcl/marshmallow.tcl]} {
   source  /home/cfa/tcl/marshmallow.tcl
}

atomselect macro dppc_head "resname DPPC and name C1 HA HB C11 H11A H11B C12 H12A H12B C13 H13A H13B H13C C14 H14A H14B H14C C15 H15A H15B H15C P O11 O12 O13 O14 N"
atomselect macro dppc_tail "resname DPPC and not name C1 HA HB C11 H11A H11B C12 H12A H12B C13 H13A H13B H13C C14 H14A H14B H14C C15 H15A H15B H15C P O11 O12 O13 O14 N"

if {[file exists /home/cfa/tcl/aln.tcl]} {
   source  /home/cfa/tcl/aln.tcl
}
if {[file exists /home/cfa/tcl/ss.tcl]} {
   source  /home/cfa/tcl/ss.tcl
}
if {[file exists /home/cfa/tcl/tmove.tcl]} {
   source /home/cfa/tcl/tmove.tcl
}

proc vzv { molid } {
    set eye_vector [vectrans \
			[measure inverse [lindex [molinfo $molid get rotate_matrix] 0]] \
			{0 0 1}
		   ]
    return $eye_vector
}

proc vxv { molid } {
    set eye_vector [vectrans \
			[measure inverse [lindex [molinfo $molid get rotate_matrix] 0]] \
			{1 0 0}
		   ]
    return $eye_vector
}

proc vyv { molid } {
    set eye_vector [vectrans \
			[measure inverse [lindex [molinfo $molid get rotate_matrix] 0]] \
			{0 1 0}
		   ]
    return $eye_vector
}

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

lappend auto_path /home/cfa/tcl/la1.0
lappend auto_path /home/cfa/tcl/orient
atomselect macro glycan {resname NAG MAN BMA FUC GAL BGNA AMAN BMAN AFUC BGAL} 
menu main on

