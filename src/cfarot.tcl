# tcl library for bond rotation
# cameron f abrams
# cfa22@drexel.edu
# 2018-2020

set CFAROT_VERSION 0.10

if {![info exists PSFGEN_BASEDIR]} {
    if {[info exists env(PSFGEN_BASEDIR)]} {
        set PSFGEN_BASEDIR $env(PSFGEN_BASEDIR)
    } else {
        set PSFGEN_BASEDIR $env(HOME)/research/psfgen
    }
}
load ${PSFGEN_BASEDIR}/lib/bondstruct.so cfa_bondstruct

###########################################################
# Primitive data handling procedures
###########################################################

# ListToArray: allocates a new array and assigns it values
# from the Tcl list $l; returns pointer to new array.  List
# elements are treated as double-precision by default.
#
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

# ListToArray_Data:  Assigns elements of an existing 
# array pointed to by $a from the tcl list $l
#
proc ListToArray_Data { a l } {
    set i 0
    foreach item $l {
	   array_setitem $a $i $item
	   incr i 1
    }
    return $a
}

# intListToArray: like listToArray except elements of
# $l are treated as integers
#
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

# intListToArray_Data: list listToArray_Data except elements
# of $l are treated as integers
#
proc intListToArray_Data {a l} {
   set length [llength $l]
   set i 0
   foreach item $l {
        arrayint_setitem $a $i $item
        incr i 1
    }
  return $a
}

# ArrayToList:  Creates a new tcl list of length $n and
# assigns it elements of existing array $a
#
proc ArrayToList {a n} {
    set length $n
    set l {}
    for {set i 0} { $i < $length} {incr i} {
	lappend l [array_getitem $a $i]
    }
    return $l
}

# intArrayToList: like ArrayToList but treats
# elements of $a as integers
#
proc intArrayToList {a n} {
    set length $n
    set l {}
    for {set i 0} { $i < $length} {incr i} {
	lappend l [arrayint_getitem $a $i]
    }
    return $l
}

proc ring_nonrotatables { ri ringsize bs } {
    if { [llength $ri] == 0 } {
        return 0
    }
    set bondsset 0
    for { set r 0 } { $r < [llength $ri] } { incr r $ringsize } {
        for {set j 0} {$j<$ringsize} {incr j} {
            set i [lindex $ri [expr $r + $j]]
            set k [lindex $ri [expr ($i+1)%$ringsize]]
            bondstruct_setbond_rotatable $bs $i $k 0
            bondstruct_setbond_rotatable $bs $k $i 0
            set bondsset [expr $bondsset + 2]
        }
    }
    return $bondsset
}

proc make_bondstruct { molid sel falist } {
    # get list of atom indices and bondlist
    puts "BONDSTRUCT) Generating bondlists for sel with [$sel num] atoms..."
    set il [$sel get index]
    set bl [$sel getbonds]
    # get number of atoms
    set na [llength $il]
    # first, count the total number of bonds in the bondlist,
    # excluding bonds to atoms outside the sel
    puts "BONDSTRUCT) Pruning bondlists..."
    flush stdout
    set bondcount 0
    for { set i 0 } { $i < $na } { incr i } {
        set a [lindex $il $i]
        set ibl [lindex $bl $i]
#        puts "$a : $ibl"
        set bb {}
        foreach pp $ibl {
            if { [lsearch $il $pp] != -1 } {
                lappend bb $pp
            }
        }
        lset bl $i $bb
        incr bondcount [llength $bb]
    }
    
    puts "BONDSTRUCT) Importing $bondcount bonds into bondstruct..."
    flush stdout
    # set up an empty bondstruct and populate it atomwise
    set ia [intListToArray $il]
    set bs [new_bondstruct $ia [llength $il] $bondcount]
#    puts "created bs with $bondcount total bonds"
    for { set i 0 } { $i < $na } { incr i } {
        set a [lindex $il $i]
        set ibl [lindex $bl $i]
        set ta [intListToArray $ibl]
        bondstruct_importbonds $bs $a $ta [llength $ibl]
#        puts "imported bonds from $a"
    }

    # any bond in a 5- or 6-membered ring, or is a peptide bond, is tagged as non-rotatable
    # any bond for which either member has as its sole heavy-atom ligand the *other* member
    # is not rotatable
    set r5 [atomselect $molid "ringsize 5 from ([$sel text])"]
    set r5i [$r5 get index]
    set n5nr [ring_nonrotatables $r5i 5 $bs]
    puts "BONDSTRUCT) Disabled rotation of $n5nr bonds in [expr [llength $r5i]/5] 5-membered rings"
    set r6 [atomselect $molid "ringsize 6 from ([$sel text])"]
    set r6i [$r6 get index]
    set n6nr [ring_nonrotatables $r6i 6 $bs]
    puts "BONDSTRUCT) Disabled rotation of $n6nr bonds in [expr [llength $r6i]/6] 6-membered rings"
    set c [atomselect $molid "protein and name C and ([$sel text])"]
    set ci [$c get index]
    set n [atomselect $molid "protein and name N and ([$sel text])"]
    set ni [$n get index]
    set nn {}
    set nnr 0
    foreach cc $ci {
        set cci [lsearch $il $cc]
        set cbl [lindex $bl $cci]
        foreach cbln $cbl {
            set nn [lsearch $ni $cbln]
            if { $nn != -1 } {
                bondstruct_setbond_rotatable $bs $cc $nn 0
                bondstruct_setbond_rotatable $bs $nn $cc 0
                set nnr [expr $nnr + 2]
                continue
            }
        }
    }
    puts "BONDSTRUCT) Disabled rotation of $nnr peptide bonds"
    flush stdout
    # any atom index which has only one bondlist member
    set nnr 0
    foreach a $il abl $bl {
        if { [llength $abl] == 1 } {
            bondstruct_setbond_rotatable $bs $a [lindex $abl 0] 0
            bondstruct_setbond_rotatable $bs [lindex $abl 0] $a 0
            set nnr [expr $nnr + 2]
        }
    }
    puts "BONDSTRUCT) Disabled rotation of $nnr single-ligand bonds"
    puts "BONDSTRUCT) Making right-sides..."
#    puts "mapping rotatables..."
#    puts "making right-sides..."
    # make the right-side lists for each bond
    bondstruct_makerightsides $bs
#    puts "printing..."
#    bondstruct_print $bs
    puts "BONDSTRUCT) Deactivating bonds with fixed atoms in their right-sides..."
    foreach fa $falist {
        bondstruct_deactivate_by_fixed $bs $fa
    }
    # generate the count of rotatable bonds and the map to the bond array
    bondstruct_maprotatables $bs
    puts "BONDSTRUCT) $bondcount total bonds, [bondstruct_getnrb $bs] rotatables"
    #puts "make_bondstruct returns"
    return $bs
}
proc bondstruct_getbond { bs i } {
   return [intArrayToList [bondstruct_getbondpointer $bs $i] 2]
}

proc bondrot_by_index { bs molid b deg } {
   set pair [bondstruct_getbond $bs $b]
   set l [lindex $pair 0]
   set r [lindex $pair 1]
   set alist [intArrayToList [bondstruct_getrightside_pointer $bs $b] [bondstruct_getrightside_count $bs $b]]
#   puts "#### bondrot_by_index attempting $b: $l - $r by $deg deg"
 #  puts "####    rightside: $alist"
   set ls [atomselect $molid "index $l"]
   set rs [atomselect $molid "index $r"]
   set lr [lindex [$ls get {x y z}] 0]
   set rr [lindex [$rs get {x y z}] 0]
   
   set asel [atomselect $molid "index $alist"]
   
#   puts "rotating [$asel num] atoms around $l - $r by $deg degrees"
   set tmat [trans bond $lr $rr $deg degrees]
   $asel move $tmat
   $ls delete
   $rs delete
   $asel delete 
}

proc bondrot_by_atomindicies { bs molid l r deg } {
   set b  [bondstruct_getbondindex $bs $l $r]
   bondrot_by_index $bs $molid $b $deg
}

# Normal IC table entry:
# I
#  \
#   \
#    J----K
#          \
#           \
#            L
# values (Rij),(Tijk),(Pijkl),(Tjkl),(Rkl)
#
# Improper type of IC table entry
# I     L
#  \   /
#   \ /
#   *K
#    |
#    |
#    J
# values (Rik),(Tikj),(Pijkl),T(jkl),(Rkl)
proc ICs_from_bondlist { molid sel } {
    set bl [$sel getbonds]
    set nm [$sel get name]
    set ai [$sel get index]
    for { set i 0 } { $i < [llength $ai] } { incr i } {
        set ilook([lindex $ai $i]) $i
        set alook([lindex $nm $i]) $i
    }
    set hs {}
    for { set j 0 } { $j < [llength $ai] } { incr j } {
        set jn [lindex $nm $j]
        set jnnn [split $jn {}]
        if { [lindex $jnnn 0] != "H" } {
            set J [lindex $ai $j]
            foreach K [lindex $bl $j] {
                # bond J-K
                set k $ilook($K)
                set kn [lindex $nm $k]
                set knnn [split $kn {}]
                if { [lindex $knnn 0 ] != "H"} {
                    set RJK [measure bond [list $J $K]]
                    foreach I [lindex $bl $j] {
                        if { $I != $K } {
                            set i $ilook($I)
                            set in [lindex $nm $i]
                            set innn [split $in {}]
                            if { ([lindex $innn 0] != "H") || ([lindex $innn 0] == "H" && [lsearch $hs $in] == -1 ) } {
                                if { [lindex $innn 0] == "H" && [lsearch $hs $in] == -1 } { lappend hs $in }
                                set RIJ [measure bond [list $I $J]]
                                set TIJK [measure angle [list $I $J $K]]
                                foreach L [lindex $bl $k] {
                                    if { $L != $J }  {
                                        set l $ilook($L)
                                        set ln [lindex $nm $l]
                                        set lnnn [split $ln {}]
                                        if { ([lindex $lnnn 0] != "H") || ([lindex $lnnn 0] == "H" && [lsearch $hs $ln] == -1 ) } {
                                            if { [lindex $lnnn 0] == "H" && [lsearch $hs $ln] == -1 } { lappend hs $ln }
                                            # I-J--K-L is a dihedral and either I or L is a yet-to-be seen H
                                            set RKL [measure bond [list $K $L]]
                                            set TJKL [measure angle [list $J $K $L]]
                                            set PIJKL [measure dihed [list $I $J $K $L]]
                                            puts -nonewline "IC [lindex $nm $i] [lindex $nm $j] [lindex $nm $k] [lindex $nm $l] "
                                            puts "[format %.4f $RIJ] [format %.4f $TIJK] [format %.4f $PIJKL] [format %.4f $TJKL] [format %.4f $RKL]"
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
