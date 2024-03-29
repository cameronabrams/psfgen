# 6VSB -- Soluble, stabilized SARS-Cov-2 S trimeric spike

## Files

This directory contains five files:
1. mkpsf_6vsb.tcl -- VMD/psfgen script that creates the first vacuum psf/pdb pair.  It uses a monte-carlo-based loop model-builder to build in the missing residues.
2. my_6vsb_vac.namd -- NAMD configuration file used to relax the "guessed" coordinates resulting from step 1.
3. my_6vsb_solv.tcl -- VMD script that uses solvate and autoionize to generate a neutralized, solvated MD system using the coordinates from step 2 as input.
4. my_6vsb_solv.namd -- NAMD configuration file that performs a minimization and short MD of the raw solvated system 
5. my_6vsb_solv_stageN.namd -- NAMD configuration file template for perform ing a minimization and series of short MD simulations of the raw solvated system.  This file is used if the `-stage` flag is set, as described below, to run the final solvated MD simulation in stages to avoid patch-grid errors.

## Instructions

Make sure PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local NAMD2 executable.  For me, these are /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

```
$ mkdir 6vsb
$ cd 6vsb
$ $PSFGEN_BASEDIR/scripts/do_test.sh -pdb 6vsb [-stage] [-npe #] [-psfgen_args [-seed #] CLEAVE WT]
```

The optional `-stage` flag, if present, instructs the script to perform the solvated MD simulations in stages to avoid patch-grid errors arising from box shrinkage during volume equilibration.  The `-npe` flag allows the user to specify the number of processing cores to use in the solvated MD simulation; 8 is the default.  The optional `-psfgen_args` flag passes subsequent arguments to the mkpsf script.  The optional `-seed` flag allows the user the specify the seed for the random-number generator.

The "CLEAVE" flag cleaves each chain at the furin cleavage site (between 685 and 686) and reverts residues at the site back to arginines.  The "WT" flag undoes the two proline mutations: P986 to K, and P987 to V.

## REFERENCE

Wrapp, Daniel; Wang, Nianshuang; Corbett, Kizzmekia S; Goldsmith, Jory A; Hsieh, Ching-Lin; Abiona, Olubukola; Graham, Barney S; McLellan, Jason S.  "Cryo-EM structure of the 2019-nCoV spike in the prefusion conformation." *Science* 2020; doi 10.1126/science.abb2507.

2017-2020, Cameron F Abrams
