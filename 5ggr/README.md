# 5GGR -- Fab from Opdivo (nivolumab) in complex with PD-1

## Files

This directory contains five files for creating a basic solvated protein system:
1. mkpsf_5ggr.tcl -- VMD/psfgen script that creates the first vacuum psf/pdb pair.  It uses a monte-carlo-based loop model-builder to build in the missing residues. 
2. my_5ggr_vac.namd -- NAMD configuration file used to relax the "guessed" coordinates resulting from step 1.
3. my_5ggr_solv.tcl -- VMD script that uses solvate and autoionize to generate a neutralized, solvated MD system using the coordinates from step 2 as input.
4. my_5ggr_colvars_op.inp -- colvars input file that defines collective variables that allow for center-of-mass restraint and an orientational restraint to keep the pseudo C2v axis of the Fab along z.
5. my_5ggr_solv.namd -- NAMD configuration file that performs a minimization and short MD of the raw solvated system; uses colvars module input file from 4.

This files are used by the generic `do_test.sh` script, as described below.

## Instructions

Make sure PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local NAMD2 executable.  For me, these are /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

The script $PSFGEN_BASEDIR/scripts/do_test.sh will build both the complex asystem and optionally, the system with just the Fab and *no* PD1.  

Option 1

```
mkdir my_5ggr
cd my_5ggr
$PSFGEN_BASEDIR/scripts/do_test.sh -pdb 5ggr
```

Option 2, no PD1:

```
mkdir my_5ggr_Fab
cd my_5ggr_Fab
$PSFGEN_BASEDIR/scripts/do_test.sh -pdb -psfgen_args -nopd1
```

## Including exipients

The script `do_sucr.sh` in this directory will build a clean system with only the Fab, waters, Cl- counterions (to neutralize the Fab), and enough surose molecules to satisfy a desired concentration in the solvent.  Here's how:

```
mkdir my_5ggr_Fab_sucr
cd my_5ggr_Fab_sucr
$PSFGEN_BASEDIR/5ggr/do_sucr.sh -cs 0.2 -psfgen_args -nopd1
```

The value after the `-cs` switch is the concentration of sucrose in Molar.  The number of sucrose molecules is computed to achieve this concentration in the non-protein-occupied volume.

Note that this requires `packmol` to be in your path.  Similar tools for other exipients are forthcoming.

2018, Cameron F Abrams
