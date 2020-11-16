# tcl library for bond rotation
# cameron f abrams
# cfa22@drexel.edu
# 2018-2020

set CFAROT_VERSION 0.10

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

# molid is the molecule id
# sel is an atomselection
# rol_il is a list of atom indices that define rotatable bond "left" partners
# rot_jl is a list of atom indices that define rotatable bond "right" partners
proc make_bondstruct { molid sel rot_il rot_jl } {

   #puts "#### rot_il $rot_il"
   #puts "#### roj_jl $rot_jl"
   set il [$sel get index]
   set bl [$sel getbonds]

   set ia [intListToArray $il]
   set bs [new_bondstruct $ia [llength $il]]
   

   foreach a $il ibl $bl {
#     if { [llength $ibl] > 0 } {
        # fix: do not add a partner atom that is not internal to the selection
        set partners [list]
        foreach pp $ibl {
            if { [lsearch $il $pp] != -1 } {
                lappend partners $pp
            }
        }
        #puts "#### $a bonds with $ibl partners $partners"
        set ta [intListToArray $partners]
        bondstruct_addbonds $bs $a $ta [llength $partners]
#     }
   }
   bondstruct_makerotatablebondlist $bs [intListToArray $rot_il] [llength $rot_il] [intListToArray $rot_jl] [llength $rot_jl]
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
   puts "#### bondrot_by_index attempting $b: $l - $r by $deg deg"
   puts "####    rightside: $alist"
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
