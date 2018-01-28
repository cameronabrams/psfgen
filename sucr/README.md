# SUCR -- a single sucrose molecule extracted from 5o8l.pdb

## Files

This directory contains five files:
1. mkpsf_sucr.tcl -- this is the VMD/psfgen script that creates the first vacuum psf/pdb pair.
2. my_sucr_vac.namd -- this is a NAMD configuration file used to relax the "guessed" coordinates resulting from step 1.
3. my_sucr_solv.tcl -- this is a VMD script that uses solvate and autoionize to generate a neutralized, solvated MD system using the coordinates from step 2 as input.
4. my_sucr_solv_stageN.namd -- this is a NAMD configuration file template for a set of MD stages that equilibrate the solvated system.  Since the system is rather small, the barostat shrinks the volume is such a way that the patch grid is invalidated early on.
5. do_test.sh -- this BASH script does everything

## Instructions

Make sure PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local namd2 executable.  For me, these are ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

1. Execute do_test.sh a clean directory

> $PSFGEN_BASEDIR/sucr/do_test.sh

The resulting files sol-stage5.coor, sol-stage5.vel, and sol-stage5.xsc are the coordinates, velocities, and extended-system files of after the final equilibration stage.

2018, Cameron F Abrams
