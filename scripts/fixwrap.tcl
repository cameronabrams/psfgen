# fixes dcd's in which one or more subunits of a multisubunit complex
# are periodically shifted but others aren't
#
# cameron f abrams cfa22@drexel.edu

# check to see if any component of a displacement vector $disp
# is outside [-hbox,hbox]
proc mycheck { disp hbox } {
    set ret [list]
    for {set i 0} {$i < 3} {incr i} {
        set a [lindex $disp $i]
        set h [lindex $hbox $i]
        set mh [expr -1*$h]
        if {$a > $h} {
            lappend ret 1
        } elseif {$a < $mh} {
            lappend ret [expr -1]
        } else {
            lappend ret 0
        }
    }
    return $ret
}

# computes all inter-subunit com-to-com displacement vectors and
# checks if any are outside [-hbox,hbox].  Identifies subunit
# move most likely to correct any violations, and returns 
# this info in a structured list to caller.  usels is a list of
# previously defined atomselections, one per subunit.
proc allchecks { molid frame usels } {
    set box [molinfo $molid get {a b c}]
    set hbox [vecscale 0.5 $box]
    set cntrs [list]
    set propmoves [list]
    foreach s $usels {
        $s frame $frame
        lappend cntrs [measure center $s]
    }
    for {set si 0} {$si < [llength $cntrs]} {incr si} {
        for {set sj 0} {$sj < [llength $cntrs]} {incr sj} {
            set thisijcheck [mycheck [vecsub [lindex $cntrs $si] [lindex $cntrs $sj]] $hbox]
            set isamove 0
            foreach w {0 1 2} {
                set diractive [lindex $thisijcheck $w]
                if { $diractive != 0 } {
                    set isamove 1
                }
            }
            if { $isamove == 1 } {
                lappend propmoves [list $sj $thisijcheck]
            }
        }
    }

    if {[llength $propmoves] > 0} {
        set counters {}
        foreach item $propmoves {
            dict incr counters $item
        }
        set maxcount 0
        set maxmove {0 0 0}
        dict for {item count} $counters {
            if { $count > $maxcount } {
                set maxcount $count
                set maxmove $item
            }
        }
        return $maxmove
    }
    return "NONE"
}

set PSF ""
set DCD ""
set SHIFTEDDCD ""

for { set i 0 } { $i < [llength $argv] } { incr i } {
    if { [lindex $argv $i] == "-psf" } {
       incr i
       set PSF [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-dcd" } {
       incr i
       set DCD [lindex $argv $i]
    }
     if { [lindex $argv $i] == "-shifteddcd" } {
       incr i
       set SHIFTEDDCD [lindex $argv $i]
    }
   if { [lindex $argv $i] == "-umacs"} {
       incr i
       set UMACS [lindex $argv $i]
    }
}

# all subunits must be defined as atomselect macros in a separate input
# file specified by the "-umacs" argument.  This file should contain the
# macro definition and a list "umacs" that just contains the list of 
# macro names.
if { [file exists $UMACS] } {
    source $UMACS
    puts "User-supplied macros in $UMACS: $umacs"
} else {
    puts "Error: No subunit definitions provided; $UMACS not found."
    exit
}

mol new $PSF
mol addfile $DCD waitfor all
set nf [molinfo top get numframes]
puts "$DCD: $nf frames"

set usels [list]
foreach um $umacs {
    lappend usels [atomselect top $um]
}

for {set f 0} { $f < $nf} {incr f} {
    animate goto $f
    set box [molinfo top get {a b c}]
    set dothismove [allchecks top $f $usels]
    if {$dothismove != "NONE"} {
        set sb [lindex $dothismove 0]
        set by [lindex $dothismove 1]
        puts "Frame $f: moving subunit [lindex $umacs $sb] by [vecmul $by $box]"
        [lindex $usels $sb] moveby [vecmul $by $box]
    } else {
        puts "Frame $f: no move suggested."
    }
}

puts "Writing subunit-shifted trajectory $SHIFTEDDCD"
animate write dcd $SHIFTEDDCD beg 0 end [expr $nf - 1] 0
exit
