# 2JIU -- Human EGFR kinase, T790M mutant, ATP-bound

## Files

This directory contains six files:
1. mkpsf_2jiu.tcl -- VMD/psfgen script that creates the first vacuum psf/pdb pair.  It uses a monte-carlo-based loop model-builder to build in the missing residues. It also replaces the ATP-like ligand with MgATP using the ADP-bound pose of 2GS6 and the Mg2+ pose relative to the last two phosphates of ATP in 5UDS.
2. my_2jiu_vac.namd -- NAMD configuration file used to relax the guessed coordinates resulting from step 1.
3. my_2jiu_solv.tcl -- VMD script that uses solvate and autoionize to generate a neutralized, solvated MD system using the coordinates from step 2 as input.
4. my_2jiu_colvars_op.inp -- colvars input file that defines collective variables that allow for center-of-mass restraint and an orientational restraint to keep the C3v axis along z.
5. my_2jiu_solv.namd -- NAMD configuration file that performs a minimization and short MD of the raw solvated system; uses colvars module input file from 4.
6. do_test.sh -- a BASH script that builds the solvated system beginning with an empty directory.

## Instructions

Make sure PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local NAMD2 executable.  For me, these are /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

1. Make a clean directory, cd into it, and issue the script command.  The files my_2jiu_i.psf and sol.coor, sol.vel, sol.xsc are the resulting system input files:

```
$ mkdir my_2jiu
$ cd my_2jiu
$ $PSFGEN_BASEDIR/2jiu/do_test.sh
``` 

2018, Cameron F Abrams, cfa22@drexel.edu
