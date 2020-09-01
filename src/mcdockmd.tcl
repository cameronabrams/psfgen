# mcdockmd
#
# Use Monte Carlo sampling of relative distance, orientation, and MD trajectory frames
# to find low-energy complexes
# 
# Modified version of script written originally by Cameron Abrams.
# The energy evaluation is being carried out in a separate C-code here
# SWIG, a Tcl-C interface building tool is used to generate `energy' package
# (c) June 2008 Harish Vashisth
# 2009, cfa
# 2017, cfa

# check for base directory name environment variable;
# if not set, use default
if {[info exists env(PSFGEN_BASEDIR)]} {
    set PSFGEN_BASEDIR $env(PSFGEN_BASEDIR)
} else {
    set HOME $env(HOME)
    set PSFGEN_BASEDIR ${HOME}/psfgen
}

# will die if CFACV_BASEDIR is invalid
source ${PSFGEN_BASEDIR}/tcl/cfa_docklib.tcl

proc main { argv } { 
    global M_PI
    global dataspace
    global HOME

    set seed 290847
    set par [list ${HOME}/charmm/toppar/par_all36_prot.prm]
    set ligdcd {"lig.psf" "ligMD.dcd" 0}
    set targdcd {"targ.psf" "targMD.dcd" 0}
    set ligsel "chain A"
    set targsel "chain A"
    set lig_reference_frame 0
    set targ_reference_frame 0
    set targ_sample_dframe 1

    set cycles 1000
    set kT 0.4 ;# kcal/mol
    set dAngleDeg 1.0
    set rmin 0.1
    set rmax 1.0
    
    set move_weights {{LIG_SAMP 0.2} {TARG_SAMP 0.4} {EULER 1.0}}

    set outputInterval 1

    set CONFIG [lindex $argv 0]

    if {[file exists $CONFIG]} {
	source $CONFIG
    } else {
	puts "ERROR: $CONFIG not found."
	exit
    }

    set dumran1 [expr srand($seed)]

    read_charmm_nb $par typ eps sig

    set targ_id [read_trajectory $targdcd]
    set targ_numframes [molinfo $targ_id get numframes]
    puts "INFO: target trajectory has $targ_numframes frames."

    set lig_id [read_trajectory $ligdcd]
    set lig_numframes [molinfo $lig_id get numframes]
    puts "INFO: ligand trajectory has $lig_numframes frames"

    set targfrm $targ_reference_frame
    set ligfrm $lig_reference_frame

    set targ [atomselect $targ_id "all" frame $targfrm]
    set lig [atomselect $lig_id "all" frame $ligfrm]
    set targref [atomselect $targ_id "all" frame $targ_reference_frame]
    set ligref [atomselect $lig_id "all" frame $lig_reference_frame]

# generate index arrays for atom types for each fragment
    set TI1 [$targ get serial]
    set i 0
    foreach t [$targ get type] {
	set ti [lsearch $typ $t]
	#    puts "$t $ti $i"
	if { [expr $ti == -1] } {
	    puts "ERROR! target atom type $t not found in parameter file."
	    exit
	}
	lset TI1 $i $ti
	incr i
    }
    
    set TI2 [$lig get serial]
    set i 0
    foreach t [$lig get type] {
	set ti [lsearch $typ $t]
	#    puts "$t $ti $i"
	if { [expr $ti == -1] } {
	    puts "ERROR! ligand atom type $t ($i) not found in parameter file."
	    exit
	}
	lset TI2 $i $ti
	incr i
    }

    # Index arrays for epsilon and sigma
    set index1 [intListToArray $TI1]
    set index2 [intListToArray $TI2]

    # Array making for Epsilon and Sigma values
    puts "$eps"
    set epsilon [ListToArray $eps]
    set sigma   [ListToArray $sig]
    # Charges arrays for both selections
    set ch1 [ListToArray [$targ get charge]]
    set ch2 [ListToArray [$lig get charge]]
    set sizeeps [llength $eps]
    set sizesig [llength $sig]
    set sizetarg [$targ num]
    set sizelig [$lig num]

    # specify number of MC cycles
    set nc [lindex $cycles 0]

    # precompute some limits based on the move magnitudes
    set zmin [expr -$dAngleDeg]
    set zmax $dAngleDeg
    set cmin 1.0
    set cmax [expr {cos($M_PI/180.0*$dAngleDeg)}]

    # Time at the beginning of Monte-Carlo trial moves loop
    puts "INFO: Docking ($nc) started at : [clock format [clock seconds] -format {%X %D} ]"

    initialize_dataspace $targ $lig

    set targR [transidentity]
    set ligR [transidentity]

    set E0 [my_call_roughenergy $targ $lig $epsilon $sigma $index1 $index2 $ch1 $ch2 $sizetarg $sizelig]

    set E $E0
    set cyclesSinceLast 0
    set ifp [open "dock.info" "w"]
    puts $ifp "# MC docking [clock format [clock seconds] -format {%X %D} ]"
#    puts $ifp "# ligand dcd: 
    puts $ifp "#  cycle  E ligframe targframe"
    set i 0
    puts $ifp "[format %08d $i]  [format %.3f $E]   [format %5d $ligfrm]    [format %5d $targfrm]"
    flush $ifp
    for {set i 1} { $i <= $nc } {incr i} {
	set chs [expr rand()]
	set MOVE [get_move $chs $move_weights]
	if {$MOVE == "TARG_SAMP"} {
	    set oldtargfrm $targfrm
	    set x [randomint 1 $targ_sample_dframe]
	    set dir 1
	    if {[randomint 0 1]} {
		set dir -1
	    }
	    set x [expr $x * $dir]
	    set targfrm [expr {$oldtargfrm+$x}]
	    if { $targfrm >= $targ_numframes } {
		set targfrm [expr {$oldtargfrm-$x}]
	    }
	    if { $targfrm < 0 } {
		set targfrm [expr {$oldtargfrm+$x}]
	    }
	    $targ frame $targfrm
	    $targref frame $oldtargfrm
	    set this_targR [measure fit $targ $targref]
	    $targ move $this_targR
	} elseif {$MOVE == "LIG_SAMP"} {
	    set oldligfrm $ligfrm
	    set x [randomint 1 $lig_sample_dframe]
	    set dir 1
	    if {[randomint 0 1]} {
		set dir -1
	    }
	    set x [expr $x * $dir]
	    set ligfrm [expr {$oldligfrm + $x}]
	    if { $ligfrm >= $lig_numframes } {
		set ligfrm [expr {$oldligfrm-$x}]
	    }
	    if { $ligfrm < 0 } {
		set ligfrm [expr {$oldligfrm+$x}]
	    }
	    $lig frame $ligfrm
	    $ligref frame $oldligfrm
	    set this_ligR [measure fit $lig $ligref]
	    $lig move $this_ligR
	} elseif {$MOVE == "EULER"} {
	    set oldpos [$lig get {x y z}]
	    # set up Euler rotation + translation trial move transformation matrix
	    set cnt [measure center $lig]
	    set z1 [random $zmin $zmax]
	    set z2 [random $zmin $zmax]
	    set cx [random $cmin $cmax]
	    set x [expr acos($cx)]
	    set z [random 0 1]
	    if { [expr $z < 0.5 ] } {
		set x [expr -$x]
	    }
	    set x [expr {$x * 57.296}]
	    # random distance
	    set r [random $rmin $rmax] 
	    # random direction
	    set rc [random -1.0 1.0]
	    set rp [expr acos($rc)]
	    set a [random -3.14159 3.14159]
	    set dx [expr {$r*sin($rp)*cos($a)}]
	    set dy [expr {$r*sin($rp)*sin($a)}]
	    set dz [expr {$r*cos($rp)}]
	    #puts "move $i: euler_zxz(deg): $z1 $x $z2; disp: $dx $dy $dz"
	    
	    set Eulerz1 [trans center $cnt axis z $z1 deg]
	    set Eulerx  [trans center $cnt axis x $x deg]
	    set Eulerz2 [trans center $cnt axis z $z2 deg]
	    set this_ligR [transmult $Eulerz1 $Eulerx $Eulerz2]
	    
	    # rotate about the COM
	    $lig move $this_ligR
	    # translate
	    $lig moveby [list $dx $dy $dz]
	} else {
	    puts "ERROR: MOVE $MOVE is not a valid move."
	    exit
	}

	set E [my_call_roughenergy $targ $lig $epsilon $sigma $index1 $index2 $ch1 $ch2 $sizetarg $sizelig]

	set X [expr rand()]
	set arg [expr {($E0-$E)/$kT}]
	if {[expr $arg < -20]} {
	    set B 0.0
	} elseif {[expr $arg > 2.8]} {
	    set B 1.1
	} else {
	    set B [expr {exp(($E0-$E)/$kT)}]
	}

	# reject
	if {[expr {$X > $B}]} { 
	    if { $MOVE == "TARG_SAMP" } {
		set targfrm $oldtargfrm
		$targ frame $targfrm
	    } elseif { $MOVE == "LIG_SAMP" } {
		set ligfrm $oldligfrm
		$lig frame $ligfrm
	    } else {
		$lig set {x y z} $oldpos
	    }
	# accept
	} else {
	    $targ writenamdbin "targetframes/[format %08d $i].coor"
	    $lig writenamdbin "ligandframes/[format %08d $i].coor"
	    if { $MOVE == "TARG_SAMP" } {
		#puts "targfrm $oldtargfrm -> $targfrm"
	    } elseif { $MOVE == "LIG_SAMP" } {
		#puts "ligfrm $oldligfrm -> $ligfrm"
		#$lig writenamdbin "ligandframes/[format %08d $i].coor"
	    } else {
		#puts " "
	    }
	    set E0 $E
	    puts $ifp "[format %08d $i]  [format %.3f $E]   [format %5d $ligfrm]     [format %5d $targfrm]"
	    flush $ifp

	}
	incr cyclesSinceLast
    }
    close $ifp

    # Time at the end of Monte-Carlo trial moves loop
    puts "INFO: Docking ended at : [clock format [clock seconds] -format {%X %D} ]"
}

main $argv

quit

