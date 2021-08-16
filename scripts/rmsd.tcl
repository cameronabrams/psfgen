# compute rmsd vs frame number and rmsf vs residue number
# Cameron F. Abrams cfa22@drexel.edu
#
# Note:  The parameters -sel-measure and -sel-align expect atomselect
# strings as values, but you can't use quoted strings with spaces because
# of VMD's (stupid) csh-based launch script.  So here, we make the
# convention that a string value for these two arguments has dashes
# in place of spaces, and split and join commands are used to remove them.

proc .. {from to} {
    if {$from >= $to} {
        for {set i $from} {$i <= $to} {incr i}    {lappend out $i}
    } else {
        for {set i $from} {$i >= $to} {incr i -1} {lappend out $i}
    }
    return $out
}

proc main { argv } {
    set psf ""
    set sel_measure ""
    set sel_align ""
    set max_frames 0
    set argc [llength $argv]
    flush stdout
    set dcd []
    set rmsd_outfile "rmsd.dat"
    set rmsf_outfile "rmsf.dat"
    set rmsf_basis "all"
    for { set i 0 } { $i < $argc } { incr i } {
        if { [lindex $argv $i] == "-psf" } {
            incr i
            set psf [lindex $argv $i]
        }
        if { [lindex $argv $i] == "-dcd" } {
            incr i
            lappend dcd [lindex $argv $i]
        }
        if { [lindex $argv $i] == "-rmsf-basis" } {
            incr i
            set rmsf_basis [lindex $argv $i]
        }
        if { [lindex $argv $i] == "-sel-measure" } {
            incr i
            set sel_measure [join [split [lindex $argv $i] -] " "]
        }
        if { [lindex $argv $i] == "-max-frames" } {
            incr i
            set max_frames [lindex $argv $i]
        }
        if { [lindex $argv $i] == "-chain-align" } {
            incr i
            set chain_align [join [split [lindex $argv $i] -] " "]
        }
        if { [lindex $argv $i] == "-rmsd-outfile" } {
            incr i
            set rmsd_outfile [lindex $argv $i]
        }
        if { [lindex $argv $i] == "-rmsf-outfile" } {
            incr i
            set rmsf_outfile [lindex $argv $i]
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
    if { [string length $sel_measure] == 0 } {
        vmdcon -err "Specify chain selection to be measured with with -sel-measure \"<string-with-dashes-for-spaces>\""
        exit
    }
    if { [string length $sel_align] == 0 } {
        vmdcon -info "Setting alignment selection equal to measure selection."
        set sel_align $sel_measure
    }

    mol new $psf
    foreach d $dcd {
        mol addfile $d waitfor all
    }
    
    set a0 [atomselect top $sel_align]
    set m0 [atomselect top $sel_measure]
    $a0 frame 0
    $m0 frame 0
    set a [atomselect top $sel_align]
    set m [atomselect top $sel_measure]
    set fp [open $rmsd_outfile "w"]
    set frame_lim [molinfo top get numframes]
    if { $max_frames > 0 && $max_frames < $frame_lim } {
        set frame_lim $max_frames
    }
    for {set i 0} {$i < $frame_lim } {incr i} {
        $a frame $i
        $m frame $i
        $a move [measure fit $a $a0]
        puts $fp "$i [format %.4f [measure rmsd $m $m0]]"
    }
    close $fp
    puts "Created $rmsd_outfile"

    if { $rmsf_basis == "CA" } {
        set rmsf_sel_measure "$sel_measure and name CA"
    } else {
        set rmsf_sel_measure "$sel_measure"
    }
    set rmsf_m [atomselect top $rmsf_sel_measure]
    set rmsf_per_atom [measure rmsf $rmsf_m first 0 last [expr $frame_lim - 1]]
    set resid_per_atom [$rmsf_m get resid]
    set mass_per_atom [$rmsf_m get mass]
    set resid [lsort -unique $resid_per_atom]
    set nresid [llength $resid]
    set mmsf_per_residue [list]
    set mass_per_residue [list]
    set rmsf_per_residue [list] 
    for {set i 0} {$i < $nresid} {incr i} {
        lappend mmsf_per_residue 0.0
        lappend mass_per_residue 0.0
        lappend rmsf_per_residue 0.0
    }
    foreach d $rmsf_per_atom r $resid_per_atom m $mass_per_atom {
        set this_msf_by_atom [expr $d*$d]
        set i [lsearch $resid $r]
        lset mmsf_per_residue $i [expr [lindex $mmsf_per_residue $i] + $m*$this_msf_by_atom]
        lset mass_per_residue $i [expr [lindex $mass_per_residue $i] + $m]
    }
    for {set i 0} {$i < $nresid} {incr i} {
        lset rmsf_per_residue $i [expr sqrt([lindex $mmsf_per_residue $i]/[lindex $mass_per_residue $i])]
    }
    set fp [open $rmsf_outfile "w"]
    foreach r $resid f $rmsf_per_residue {
        puts $fp "$r [format %.4f $f]"
    }
    close $fp
    vmdcon -info "Created $rmsf_outfile"

}

main $argv
exit


