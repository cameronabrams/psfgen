# cameron f abrams (c) 2019
# drexel university
# chemical and biological engineering
#
# check for base directory name variable;
# if not set, use default
if {![info exists PSFGEN_BASEDIR]} {
  # see if user set an environment variable
  if {[info exists env(PSFGEN_BASEDIR)]} {
      set PSFGEN_BASEDIR $env(PSFGEN_BASEDIR)
  } else {
      set PSFGEN_BASEDIR $env(HOME)/research/psfgen
  }
}
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
set X ALA 
for { set a 0 } { $a < [llength $argv] } { incr a } {
   set arg [lindex $argv $a]
   if { $arg == "-x" } {
      incr a
      set X [lindex $argv $a]
   }
}

puts "$X"

package require psfgen
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf

mol new skel.pdb
set n2pos [cacoIn_nOut 1 A 0]

segment A {
  pdb skel.pdb
  residue 2 ${X} A
  residue 3 GLY A
}

coordpdb skel.pdb A
coord A 2 N $n2pos

guesscoord

writepsf "my_gxg.psf"
writepdb "my_gxg.pdb"

exit
