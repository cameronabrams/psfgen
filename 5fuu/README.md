# 5FUU -- Soluble, cleaved JR-FL trimer with glycans

## Files

This directory contains five files:
1. mkpsf_5fuu.tcl -- VMD/psfgen script that creates the first vacuum psf/pdb pair.  It uses a monte-carlo-based loop model-builder to build in the missing residues.  It includes all glycans present in the PDB.
2. my_5fuu_vac.namd -- NAMD configuration file used to relax the "guessed" coordinates resulting from step 1.
3. my_5fuu_solv.tcl -- VMD script that uses solvate and autoionize to generate a neutralized, solvated MD system using the coordinates from step 2 as input.
4. my_5fuu_colvars_op.inp -- colvars input file that defines collective variables that allow for center-of-mass restraint and an orientational restraint to keep the C3v axis along z.
5. my_5fuu_solv.namd -- NAMD configuration file that performs a minimization and short MD of the raw solvated system 
5. my_5fuu_solv_stageN.namd -- NAMD configuration file template for perform ing a minimization and series of short MD simulations of the raw solvated system.  This file is used if the `-stage` flag is set, as described below, to run the final solvated MD simulation in stages to avoid patch-grid errors.

## Instructions

Make sure PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local NAMD2 executable.  For me, these are /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

```
$ mkdir 5fuu
$ cd 5fuu
$ $PSFGEN_BASEDIR/scripts/do_test.sh -pdb 5fuu [-stage] [-npe #] [-psfgen_args [-seed #] [-mper-extend] [-tm-extend] [-log-dcd <filename>] [-man9 # [-man9 #...]]]
```

The optional `-stage` flag, if present, instructs the scripts to perform the solvated MD simulations in stages to avoid patch-grid errors arising from box shrinkage during volume equilibration.  The `-npe` flag allows the user to specify the number of processing cores to use in the solvated MD simulation; 8 is the default.  The optional `-psfgen_args` flag passes subsequent arguments to the mkpsf script.  The optional `-seed` flag allows the user the specify the seed for the random-number generator.  The optional `-mper-extend` flag, if present, instructs the script to model-in the MPER residues up to and including 682 as an alpha-helix.  The optional `-tm-extend` flag, if present, instructs the script to model-in the TM residues up to 709 as an alpha-helix; note that `-tm-extend` will set `-mper-extend` by default.  The `-log-dcd` switch instructs the script to log configurations during the build to DCD file named in the switch.  Each `-man9` flag, if present, allows the user to indicate which asparagine residues to model-build a man9 glycan onto; currently only residues 386 and 392 on gp120 are supported with this functionality.

2017-2018, Cameron F Abrams
