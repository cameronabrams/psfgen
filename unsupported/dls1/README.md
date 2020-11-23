# B529 -- a single DLS1 linker-1 molecule extracted from AEG (internal)

## Files

This directory contains six files:
1. dls1_raw.pdb -- PDB file exported from Avogadro and used to generate dls1.str using cgenff, with atoms renamed to conform to automatic naming from cgenff 
2. mkpsf_dsl1.tcl -- this is the VMD/psfgen script that creates the first vacuum psf/pdb pair.
2. my_dls1_vac.namd -- this is a NAMD configuration file used to relax the "guessed" coordinates resulting from step 1.
3. my_dls1_solv.tcl -- this is a VMD script that uses solvate and autoionize to generate a neutralized, solvated MD system using the coordinates from step 2 as input.
4. my_dls1_solv_stageN.namd -- this is a NAMD configuration file template for a set of MD stages that equilibrates the solvated system.  Since the system is rather small, the barostat shrinks the volume is such a way that the patch grid is invalidated early on.
5. V.gp -- gnuplot script for plotting volume vs timestep during to solvated MD
6. do_test.sh -- this BASH script does everything

## Instructions

Make sure PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local namd2 executable.  For me, these are ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

1. Execute do_test.sh a clean directory

```
$ mkdir dls1
$ cd dls1
$ $PSFGEN_BASEDIR/dls1/do_test.sh
```

The resulting files sol-stage2.coor, sol-stage2.vel, and sol-stage2.xsc are the coordinates, velocities, and extended-system files of after the final equilibration stage.  The final coordinates of just the DLS1 molecule are stored in `my_dls1_x.pdb`.

This uses the file dls1.str in $PSFGEN_BASEDIR/charmm, which was generated using the CGenFF server from a PDB of DLS1 from drawn in Avogadro.

2018, Cameron F Abrams
