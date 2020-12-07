# compute rmsd vs frame number
set setstr "protein and name CA"
set outfile "rmsd.dat"

set argc [llength $argv]
if { $argc < 2 } {
    puts "Error: You should specify a psf and at least one dcd as args"
    exit
}

set psf [lindex $argv 0]
if { ! [file exists $psf] } {
    puts "Error: PSF file $psf not found."
    exit
}
set dcd {}
for { set a 1 } { $a < [llength $argv] } { incr a } {
    set d [lindex $argv $a]
    lappend dcd $d
    if { ! [file exists $d] } {
        puts "Error: DCD file $d not found."
        exit
    }
}

mol new $psf
foreach d $dcd {
    mol addfile $d waitfor all
}

set a0 [atomselect top $selstr]
$a0 frame 0
set a [atomselect top "protein and name CA"]
set fp [open $outfile "w"]
for {set i 0} {$i < [molinfo top get numframes]} {incr i} {
    $a frame $i
    puts $fp "$i [format %.4f [measure rmsd $a $a0]]"
}
close $fp
exit


