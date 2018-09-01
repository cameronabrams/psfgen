#!/bin/bash
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

SYSNAME=5fuu-mvn-davei
vmd -dispdev text -e $PSFGEN_BASEDIR/${SYSNAME}/mkpsf.tcl > psfgen1.log
cp $PSFGEN_BASEDIR/${SYSNAME}/my_complex_vac_stage0.namd .
$CHARMRUN +p2 $NAMD2 my_complex_vac_stage0.namd > vac0.log
cp $PSFGEN_BASEDIR/${SYSNAME}/my_complex_vac_stage1.namd .
$CHARMRUN +p2 $NAMD2 my_complex_vac_stage1.namd > vac1.log
echo "Done."

