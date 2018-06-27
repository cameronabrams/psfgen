# tcl library for bond rotation
# cameron f abrams
# cfa22@drexel.edu
# 2018

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

proc make_bondstruct { molid sel } {
   set il [$sel get index]
   set bl [$sel getbonds]

   set ia [intListToArray $il]
   set bs [new_bondstruct $ia [llength $il]]
   
   foreach a $il ibl $bl {
     if { [llength $ibl] > 0 } {
        set ta [intListToArray $ibl]
        bondstruct_addbondlist $bs $a $ta [llength $ibl]
     }
   }
   
   return $bs
}

proc bondstruct_getrsel_indices { bs i j } {
   puts "$i $j"
   set resarray [bondstruct_getrl $bs $i $j]
   set tmp [intArrayToList $resarray [bondstruct_getna $bs]]
   set ret {}
   foreach i $tmp {
     if { $i != -1 } {
       lappend ret $i
     }
   }
   return $ret
}

proc my_bondrot { molid rsel i j deg } {
   set is [atomselect $molid "index $i"]
   set js [atomselect $molid "index $j"]
   set ir [lindex [$is get {x y z}] 0]
   set jr [lindex [$js get {x y z}] 0]
   $rsel move [trans bond $ir $jr $deg degrees]
}

proc gentwist { molid sel i j deg } {

   set bs [make_bondstruct $molid $sel]

#   puts "calling getrsel for $i $j"

   set rl [bondstruct_getrsel_indices $bs $i $j]
#   puts "rl $rl"

   set rsel [atomselect $molid "index $rl"]
   my_bondrot $molid $rsel $i $j $deg 
   
#   set names [$rsel get name]
#   puts "names $names"
   
#   print_bondlist $bs   

   

   free_bondstruct $bs

}
