# 5U7O -- SOSIP trimer with glycans and entry inhibitor BMS-529 bound

## Files

This directory contains five files:
1. mkpsf_5u7o.tcl -- the VMD/psfgen script that creates the first vacuum psf/pdb pair.  It uses a monte-carlo-based loop model-builder to build in missing residues.  It includes all glycans present in the PDB entry.  Symmetry replication is peformed by aligning the gp120/gp41 alpha carbons to the 4ZMJ structure and then using the BIOMT transformations in 4ZMJ (since 5U7O has none!) to generate the full trimer.
2. my_5u7o_vac.namd -- NAMD configuration file used to relax the "guessed" coordinates resulting from mkpsf_5u7o.tcl.
3. my_5u7o_solv.tcl -- VMD script that uses `solvate` and `autoionize` to generate a neutralized, solvated MD system using the coordinates from step 2 as input.  This also generates my_5u7o_colvars_op.inp, a colvars input file that defines collective variables that allow for center-of-mass restraint and an orientational restraint to keep the C3v axis along z.  (Colvars restraints are NOT recommended for production MD.)
4. my_5u7o_solv.namd -- NAMD configuration file that performs a minimization and short MD of the raw solvated system; uses colvars module input file from 4.
5. do_test.sh -- a bash script that does everything

## Instructions

If you have cloned this repository, then make sure PSFGEN_BASEDIR resolves to the root directory of your local copy.  If you did not
clone the repository, you will have to figure this part out on your own.  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local NAMD2 executable.  For me, these are /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

To create the system, just invoke the do_test script in a clean directory:

```
$ mkdir 5u7o
$ cd 5u7o
$ $PSFGEN_BASEDIR/5u7o/do_test.sh
```

2018, Cameron F Abrams
