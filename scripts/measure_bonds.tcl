# measure_bonds.tcl -- 
# reads lines from input file, each is CHAIN, RES1, RES2
# measures distance between C on RES1 and N on RES2
# appends that distance to the lines
# writes all lines back out to file

# first argument is PSF, second is PDB or COOR, third is name of input file
vmdcon -info "measure_bonds.tcl: args: $argv"
if { [llength $argv] != 4 } {
    vmdcon -error "measure_bonds.tcl expects four positional arguments:\n        psf coor infile-to-be-modified fixed-pdb"
    exit 1
}
set psf [lindex $argv 0]
set coor [lindex $argv 1]
set infile [lindex $argv 2]
set fixed [lindex $argv 3]

mol new $psf
mol addfile $coor
set a [atomselect top "all"]
$a set occupancy 0

set fp [open $infile "r"]
set lines [split [read $fp] \n]
close $fp
set RES1 {}
set RES2 {}
set CH {}
foreach l $lines {
    if { [llength $l] > 0 } {
       lappend CH [lindex $l 0]
       lappend RES1 [lindex $l 1]
       lappend RES2 [lindex $l 2]
    }
}
puts "CH $CH"
puts "RES1 $RES1"
puts "RES2 $RES2"

set bl {}
set CC {}
set NN {}
foreach c $CH i $RES1 j $RES2 {
    set ci [string index $i end]
    if {[string is alpha $ci]} {
        set resid [string range $i 0 end-1]
        set theC [atomselect top "chain $c and resid $resid and insertion $ci and name C"]
    } else {
        set theC [atomselect top "chain $c and resid $i and name C"]
    }
    set ni [string index $j end]
    if {[string is alpha $ni]} {
        set resid [string range $j 0 end-1]
        set theN [atomselect top "chain $c and resid $resid and insertion $ni and name N"]
    } else {
        set theN [atomselect top "chain $c and resid $j and name N"]
    }
    puts "resid $i on chain $c has [$theC num] Cs with serial [$theC get serial]"
    puts "resid $j on chain $c has [$theN num] Ns with serial [$theN get serial]"
    set ii [$theC get index]
    set jj [$theN get index]
    $theN set occupancy 1
    lappend bl [measure bond [list $ii $jj]]
    lappend CC [$theC get serial]
    lappend NN [$theN get serial]
}
puts "CC $CC"
puts "NN $NN"

set fp [open $infile "w"]
foreach c $CH i $CC j $NN b $bl {
    puts $fp "$c $i $j [format %.4f $b]"
}
close $fp
$a writepdb $fixed
vmdcon -info "measure_bonds.tcl: $infile modified and $fixed created"
exit
