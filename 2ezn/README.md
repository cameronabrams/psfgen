# 2ezn -- Cyanovirin-N ("CVN") with option to make CVN-(G4S)n-H6-MPER DAVEI

## Files

This directory contains five files:
1. mkpsf_2ezn.tcl -- the VMD/psfgen script that creates the first vacuum psf/pdb pair.  It uses a monte-carlo-based loop model-builder to build in missing residues.  It includes all glycans present in the PDB entry.  Symmetry replication is peformed by aligning the gp120/gp41 alpha carbons to the 4ZMJ structure and then using the BIOMT transformations in 4ZMJ (since 5U7O has none!) to generate the full trimer.
2. my_2ezn_vac.namd -- NAMD configuration file used to relax the "guessed" coordinates resulting from mkpsf_2ezn.tcl.
3. my_2ezn_solv.tcl -- VMD script that uses `solvate` and `autoionize` to generate a neutralized, solvated MD system using the coordinates from step 2 as input.  This also generates my_2ezn_colvars_op.inp, a colvars input file that defines collective variables that allow for center-of-mass restraint and an orientational restraint to keep the C3v axis along z.  (Colvars restraints are NOT recommended for production MD.)
4. my_2ezn_solv_stageN.namd -- NAMD configuration file template for performs a minimization and short NpT MD of the raw solvated system.
5. do_test.sh -- a bash script that does everything

## Instructions

If you have cloned this repository, then make sure PSFGEN_BASEDIR resolves to the root directory of your local copy.  If you did not
clone the repository, you will have to figure this part out on your own.  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local NAMD2 executable.  For me, these are /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

To create the system, just invoke the do_test script in a clean directory:

```
$ mkdir 2ezn
$ cd 2ezn
$ $PSFGEN_BASEDIR/2ezn/do_test.sh
```

## DAVEI Option

To generate a DAVEI molecule (a fusion of CVN, N repeats of (Gly4Ser), a His-tag, and the 27-residue HIV-1 MPER segment), do this:

```
# mkdir my_davei
# cd my_davei
$ $PSFGEN_BASEDIR/2ezn/do_test.sh -psfgen_args -davei <N>
```

where `<N>` is replaced by the number of repeats of the (Gly4Ser) linker unit you want.

2018, Cameron F Abrams
