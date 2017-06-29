# 1F7A -- dimeric HIV-1 protease with bound substrate

## Files

This directory contains five files:
1. mkpsf_1f7a.tcl -- this is the VMD/psfgen script that creates the first vacuum psf/pdb pair.  It uses a monte-carlo-based loop model-builder to build in the missing residues. It handles the two protonated aspartates and the two down-puckered prolines.
2. my_1f7a_vac.namd -- this is a NAMD configuration file used to relax the "guessed" coordinates resulting from step 1.
3. my_1f7a_solv.tcl -- this is a VMD script that uses solvate and autoionize to generate a neutralized, solvated MD system using the coordinates from step 2 as input.
4. my_1f7a_colvars_op.inp -- a colvars input file that defines collective variables that allow for center-of-mass restraint and an orientational restraint.
5. my_1f7a_solv.namd -- this is a NAMD configuration file that performs a minimization and short MD of the raw solvated system; uses colvars module input file from 4.

This file do_test.sh in the scripts/ directory is a Bash script that performs the sequence of commands detailed in the instructions below.

## Instructions

Make sure PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local namd2 executable.  For me, these are ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

1. Download 1f7a.pdb to a clean directory

> wget http://www.rcsb.org/pdb/files/1f7a.pdb

2. Use VMD in text mode to generate the psf/pdb

> vmd -dispdev text -e $PSFGEN_BASEDIR/1f7a/mkpsf_1f7a.tcl

3. Run NAMD to relax bonds and guessed-in atoms

> ln -s $PSFGEN_BASEDIR/1f7a/my_1f7a_vac.namd .

> $CHARMRUN +p1 $NAMD2 my_1f7a_vac.namd > vac.log

4. Use VMD to solvate and neutralize the output of step 3

> vmd -dispdev text -e $PSFGEN_BASEDIR/1f7a/my_1f7a_solv.tcl

5. Run NAMD to minimize and shake out the solvated system

> ln -s $PSFGEN_BASEDIR/1f7a/my_1f7a_trimer_solv.namd .

> $CHARMRUN +p16 $NAMD2 my_1f7a_solv.namd > solv.log

2017, Cameron F Abrams
