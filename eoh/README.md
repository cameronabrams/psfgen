# EOH -- a single ethanol molecule extracted from 3tod.pdb

## Files

This directory contains seven files:
1. mkpsf_mtl.tcl -- this is the VMD/psfgen script that creates the first vacuum psf/pdb pair.
2. my_mtl_vac.namd -- this is a NAMD configuration file used to relax the "guessed" coordinates resulting from step 1.
3. my_mtl_solv.tcl -- this is a VMD script that uses solvate and autoionize to generate a neutralized, solvated MD system using the coordinates from step 2 as input.
4. my_mtl_solv_stageN.namd -- this is a NAMD configuration file template for a set of MD stages that equilibrate the solvated system.  Since the system is rather small, the barostat shrinks the volume is such a way that the patch grid is invalidated early on.
5. extractpdb.tcl -- a VMD script that generates a packmol-ready PDB for mannitol with segname Q, chain Q.
6. V.gp -- gnuplot script for plotting volume vs timestep during to solvated MD
7. do_test.sh -- this BASH script does everything

## Instructions

Make sure PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local namd2 executable.  For me, these are ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

1. Execute do_test.sh a clean directory

```
$ mkdir eoh
$ cd eoh
$ $PSFGEN_BASEDIR/eoh/do_test.sh
```

The resulting files sol-stage2.coor, sol-stage2.vel, and sol-stage2.xsc are the coordinates, velocities, and extended-system files of after the final equilibration stage.  The final coordinates of just the ethanol molecule are stored in `my_eoh_q.pdb`.

2019, Cameron F Abrams
