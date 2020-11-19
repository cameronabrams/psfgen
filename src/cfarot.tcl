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
    for { set r 0 } { $r < [llength $ri] } { incr r $ringsize } {
        for {set j 0} {$j<$ringsize} {incr j} {
            set i [lindex $ri [expr $r + $j]]
            set k [lindex $ri [expr ($i+1)%$ringsize]]
            bondstruct_setbond_rotatable $bs $i $k 0
            bondstruct_setbond_rotatable $bs $k $i 0
        }
    }
}
# a and b are atom indices of a bond
# ri is an ordered list of ring-atom indices; every block of $ringsize
# entries is one unique ring
# this procedure returns 1 if both a and b are located in any one 
# ring
proc bond_in_ring { a b ri ringsize } {
#    puts "$a $b $ri $ringsize"
    if { [llength $ri] == 0 } {
        return 0
    }
    for { set i 0 } { $i < [llength $ri] } { incr i $ringsize } {
        set this_ring {}
        for {set j 0} {$j<$ringsize} {incr j} {
            lappend this_ring [lindex $ri [expr $i + $j]]
        }
#       puts "searching ($this_ring) for $a-$b"
#        flush stdout
        if { [lsearch $this_ring $a] != -1 && [lsearch $this_ring $b] != -1 } {
            return 1
        }
    }
    return 0
}

proc bond_is_peptide { a b ci ni } {
    if { [llength $ci] == 0 } {
        return 0
    }
    if { [lsearch $ci $a] !=-1 && [lsearch $ni $b] != -1 } {
        return 1
    }
    if { [lsearch $ci $b] !=-1 && [lsearch $ni $a] != -1 } {
        return 1
    }
    return 0
}

proc make_bondstruct { molid sel } {
    # get list of atom indices and bondlist
    set il [$sel get index]
    set bl [$sel getbonds]
    # get number of atoms
    set na [llength $il]
    # first, count the total number of bonds in the bondlist,
    # excluding bonds to atoms outside the sel
    puts "BONDSTRUCT) Pruning from VMD sel getbonds..."
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
    
    puts "BONDSTRUCT) Importing into bondstruct..."
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
    ring_nonrotatables $r5i 5 $bs
    puts "BONDSTRUCT) Considering [expr [llength $r5i]/5] 5-membered rings"
    set r6 [atomselect $molid "ringsize 6 from ([$sel text])"]
    set r6i [$r6 get index]
    ring_nonrotatables $r6i 6 $bs
    puts "BONDSTRUCT) Considering [expr [llength $r6i]/6] 6-membered rings"
    set c [atomselect $molid "protein and name C and ([$sel text])"]
    set ci [$c get index]
    set n [atomselect $molid "protein and name N and ([$sel text])"]
    set ni [$n get index]
    puts "BONDSTRUCT) Considering [llength $ci] peptide bonds"
    puts "BONDSTRUCT) Labeling rotatables..."
    flush stdout
    foreach a $il ibl $bl {
        foreach b $ibl {
#            puts "Considering $a-$b"
#            flush stdout
            set rotatable 1
            set ispeptidebond [bond_is_peptide $a $b $ni $ci]
#            puts " -> 5 $in5ring 6 $in6ring p $ispeptidebond"
            if { $ispeptidebond == 1 } {
                set rotatable 0
            }
            if { [llength $ibl] == 1 } {
                set rotatable 0
            }
            bondstruct_setbond_rotatable $bs $a $b $rotatable
    #        puts "-> $a $b $rotatable"
        }
    }
    puts "BONDSTRUCT) Making right-sides..."
#    puts "mapping rotatables..."
    # generate the count of rotatable bonds and the map to the bond array
    bondstruct_maprotatables $bs
#    puts "making right-sides..."
    # make the right-side lists for each bond
    bondstruct_makerightsides $bs
#    puts "printing..."
#    bondstruct_print $bs

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
