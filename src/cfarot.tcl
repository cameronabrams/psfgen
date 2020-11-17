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

# a and b are atom indices of a bond
# ri is an ordered list of ring-atom indices; every block of $ringsize
# entries is one unique ring
# this procedure returns 1 if both a and b are located in any one 
# ring
proc bond_in_ring { a b ri ringsize } {
    for { set i 0 } { $i < [llength $ri] } { incr i $ringsize } {
        set this_ring {}
        for {set j 0} {$j<$ringsize} {incr $j} {
            lappend this_ring [expr $i + $j]
        }
        if { [lsearch $this_ring $a] != -1 && [lsearch $this_ring $b] != -1] } {
            return 1
        }
    }
    return 0
}

proc bond_is_peptide { a b ci ni } {
    if { [lsearch $ci $a] !=-1 and [lsearch $ni $b] != -1 } {
        return 1
    }
    if { [lsearch $ci $b] !=-1 and [lsearch $ni $a] != -1 } {
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
    set bondcount 0
    for { set i 0 } { $i < $na } { incr i } {
        set a [lindex $il $i]
        set ibl [lindex $bl $i]
        puts "$a : $ibl"
        set bb {}
        foreach pp $ibl {
            if { [lsearch $il $pp] != -1 } {
                lappend partners $pp
            }
        }
        lset bl $i $partners
        incr bondcount [llength $partners]
    }
    
    # set up an empty bondstruct and populate it atomwise
    set ia [intListToArray $il]
    set bs [new_bondstruct $ia [llength $il] $bondcount]
    for { set i 0 } { $i < $na } { incr i } {
        set a [lindex $il $i]
        set ibl [lindex $bl $i]
        set ta [intListToArray $ibl]
        bondstruct_importbonds $bs $a $ta [llength $ibl]
    }

    # any bond in a 5- or 6-membered ring, or is a peptide bond, is tagged as non-rotatable
    # any bond for which either member has as its sole heavy-atom ligand the *other* member
    # is not rotatable
    set r5 [atomselect $mol "ringsize 5 from ([$sel str])"]
    set r5i [$r5 get index]
    set r6 [atomselect $mol "ringsize 6 from ([$sel str])"]
    set r6i [$r6 get index]
    set c [atomselect $molid "protein and name C and ([$sel str])"]
    set ci [$c get index]
    set n [atomselect $molid "protein and name N and ([$sel str])"]
    set ni [$n get index]
    foreach a $il ibl $bl {
        foreach b $ibl {
            set rotatable 1
            set in5ring [bond_in_ring $a $b $r5i 5]
            set in6ring [bond_in_ring $a $b $r6i 6]
            set ispeptidebond [bond_is_peptide $a $b $ni $ci]
            if { $in5ring == 1 || $in6ring == 1 || $ispeptidebond == 1 } {
                set rotatable 0
            }
            if { [llength $ibl] == 1 } {
                set rotatable 0
            }
            if { $rotatable == 0 } {
                bondstruct_setbond_rotatable $bs $a $b 0
            } else {
                bondstruct_setbond_rotatable $bs $a $b 1
            }
        }
    }
    # generate the count of rotatable bonds and the map to the bond array
    bondstruct_maprotatables $bs
    # make the right-side lists for each bond
    bondstruct_makerightsides $bs
    bondstruct_print $bs

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
   #puts "#### bondrot_by_index attempting $b: $l - $r by $deg deg"
   #puts "####    rightside: $alist"
   set ls [atomselect $molid "index $l"]
   set rs [atomselect $molid "index $r"]
   set lr [lindex [$ls get {x y z}] 0]
   set rr [lindex [$rs get {x y z}] 0]
   
   set asel [atomselect $molid "index $alist"]
   
   #puts "rotating [$asel num] atoms around $l - $r by $deg degrees"
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
