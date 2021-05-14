# compute rmsd vs frame number and rmsf vs residue number
# Cameron F. Abrams cfa22@drexel.edu
#
# Note:  The parameters -sel-measure and -sel-align expect atomselect
# strings as values, but you can't use quoted strings with spaces because
# of VMD's (stupid) csh-based launch script.  So here, we make the
# convention that a string value for these two arguments has dashes
# in place of spaces, and split and join commands are used to remove them.

proc main { argv } {
    set psf ""
    set sel_measure ""
    set sel_align ""
    set argc [llength $argv]
    puts "$argv"
    flush stdout
    set dcd []
    set rmsd_outfile "rmsd.dat"
    set rmsf_outfile "rmsf.dat"
    for { set i 0 } { $i < $argc } { incr i } {
        if { [lindex $argv $i] == "-psf" } {
            incr i
            set psf [lindex $argv $i]
        }
        if { [lindex $argv $i] == "-dcd" } {
            incr i
            lappend dcd [lindex $argv $i]
        }
        if { [lindex $argv $i] == "-sel-measure" } {
            incr i
            set sel_measure [join [split [lindex $argv $i] -] " "]
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
        puts "Error: Specify PSF with -psf"
        exit
    }
    if { [llength $dcd] == 0 } {
        puts "Error: Specify one or more DCD's with -dcd file1.dcd -dcd file2.dcd ..."
        exit
    }
    if { [string length $sel_measure] == 0 } {
        puts "Error: Specify chain selection to be measured with with -sel-measure \"<string-with-dashes-for-spaces>\""
        exit
    }
    if { [string length $sel_align] == 0 } {
        puts "Setting alignment selection equal to measure selection."
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
    for {set i 0} {$i < [molinfo top get numframes]} {incr i} {
        $a frame $i
        $m frame $i
        $a move [measure fit $a $a0]
        puts $fp "$i [format %.4f [measure rmsd $m $m0]]"
    }
    close $fp
    puts "Created $rmsd_outfile"

    set rmsf_sel_measure "$sel_measure and name CA"
    puts "$rmsf_sel_measure"
    set rmsf_m [atomselect top $rmsf_sel_measure]
    set rmsf [measure rmsf $rmsf_m]
    set fp [open $rmsf_outfile "w"]
    foreach r [$rmsf_m get resid] f $rmsf {
        puts $fp "$r [format %.4f $f]"
    }
    close $fp
    puts "Created $rmsf_outfile"

}
puts "$argv"
main $argv
exit


