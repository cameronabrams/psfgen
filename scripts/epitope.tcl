# Computes and returns a dictionary of chainA-residue occupancy frequency
# in the chainA-chainB epitope

proc update_resid_set { the_set new_set } {
    foreach r $new_set {
        set posn [lsearch $the_set $r]
        if { $posn == -1 } {
            lappend the_set $r
        }
    }
    return $the_set
}

proc update_resid_counts { counts the_set } {
    foreach r $the_set {
        set old_count 0
        if {[dict exists $counts $r]} {
            set old_count [dict get $counts $r]
        }
        set new_count [expr $old_count + 1]
        dict set counts $r $new_count
    }
    return $counts
}

proc get_epitope { chainA chainB molid } {
    set epitope []
    set counts [dict create]
    set numframes [molinfo $molid get numframes]
    set r [atomselect $molid "protein and chain $chainA and same residue as (within 3.0 of (protein and chain $chainB))"]
    for { set i 0 } { $i < $numframes } { incr i } {
        $r frame $i 
        $r update
        set this_epitope [lsort -unique [$r get resid]]
        set epitope [update_resid_set $epitope $this_epitope]
        set counts [update_resid_counts $counts $this_epitope]
        puts "-> frame ${i}..."
    }
    set occ_dict [dict create]
    set epitope [lsort $epitope]
    foreach r $epitope {
        dict set occ_dict $r [format "%.4f" [expr [dict get $counts $r]*1.0/$numframes]]
    }
    return $occ_dict
}

proc save_epitope {filename occ_dict} {
    set fp [open $filename "w"]
    dict for {r o} $occ_dict {
        puts $fp "$r $o"
    }
    close $fp
}
proc debug {} {
    set test_set {0 1 2 3 4}
    set counts [dict create]
    foreach r $test_set {
        dict set counts $r 1
    }
    set new_set { 2 3 4 5 }
    set test_test [update_resid_set $test_set $new_set]
    set counts [update_resid_counts $counts $test_set ]

    foreach r $test_set {
        vmdcon -info "resid $r count [dict get $counts $r]"
    }
}

proc main { argv } {
    set psf ""
    set chainA ""
    set chainB ""
    set max_frames 0
    set argc [llength $argv]
    flush stdout
    set dcd []
    set outfile "epitope.dat"
    for { set i 0 } { $i < $argc } { incr i } {
        if { [lindex $argv $i] == "-psf" } {
            incr i
            set psf [lindex $argv $i]
        }
        if { [lindex $argv $i] == "-dcd" } {
            incr i
            lappend dcd [lindex $argv $i]
        }
        if { [lindex $argv $i] == "-chainA" } {
            incr i
            set chainA [lindex $argv $i]
        }
        if { [lindex $argv $i] == "-chainB" } {
            incr i
            set chainB [lindex $argv $i]
        }        
        if { [lindex $argv $i] == "-o" } {
            incr i
            set outfile [lindex $argv $i]
        }
    }
    
    if { [string length $psf] == 0 } {
        vmdcon -err "Specify PSF with -psf"
        exit
    }
    if { [llength $dcd] == 0 } {
        vmdcon -err "Specify one or more DCD's with -dcd file1.dcd -dcd file2.dcd ..."
        exit
    }
    if { [string length $chainA] == 0 } {
        vmdcon -err "Specify chainA designation with -chainA \"<letter>\""
        exit
    }
    if { [string length $chainB] == 0 } {
        vmdcon -err "Specify chainB designation with -chainB \"<letter>\""
        exit
    }

    mol new $psf
    foreach d $dcd {
        mol addfile $d waitfor all
    }

    set occ [get_epitope $chainA $chainB top]
    save_epitope $outfile $occ 
}

main $argv
exit
