# 2K7W -- BAX with bound BIM SAHB peptide

## Files

This directory contains five files:
1. mkpsf_2k7w.tcl -- this is the VMD/psfgen script that creates the first vacuum psf/pdb pair.  It uses a monte-carlo-based loop model-builder to build in the missing residues. It handles the two protonated aspartates and the two down-puckered prolines.
2. my_2k7w_vac.namd -- this is a NAMD configuration file used to relax the "guessed" coordinates resulting from step 1.
3. my_2k7w_solv.tcl -- this is a VMD script that uses solvate and autoionize to generate a neutralized, solvated MD system using the coordinates from step 2 as input.
4. my_2k7w_colvars_op.inp -- a colvars input file that defines collective variables that allow for center-of-mass restraint and an orientational restraint.
5. my_2k7w_solv.namd -- this is a NAMD configuration file that performs a minimization and short MD of the raw solvated system; uses colvars module input file from 4.

The file do_test.sh in the scripts/ directory is a general Bash script that will generate a run-ready solvated system from a clean directory.

## Instructions

Make sure PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local namd2 executable.  For me, these are ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

```
mkdir my_bax
cd my_bax
${PSFGEN_BASEDIR}/scripts/do_test.sh -pdb 2k7w [-psfgen_args -bim -frm #]
```

The `-bim` flag signals to include the BIM peptide (otherwise it is not included by default).  The `-frm` flag allows you to specify which of the 20 frames in the NMR PDB file you would like to use; by default, it uses the last frame.

2019, Cameron F Abrams
