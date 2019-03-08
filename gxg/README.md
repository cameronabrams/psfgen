# GXG -- a single GXG tripeptide molecule

## Instructions

Make sure PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local namd2 executable.  For me, these are ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

1. Execute do_test.sh a clean directory

```
$ mkdir gxg
$ cd gxg
$ $PSFGEN_BASEDIR/gxg/do_test.sh -psfgen_args -x ALA
```
You can replace the 'ALA' with any three-letter amino acid residue in the CHARMM topology.

The resulting files sol-stage2.coor, sol-stage2.vel, and sol-stage2.xsc are the coordinates, velocities, and extended-system files of after the final equilibration stage.  The final coordinates of just the GXG molecule are stored in `my_gxg_q.pdb`.

2019, Cameron F Abrams
