# 5GGR -- Fab from Opdivo (nivolumab) in complex with PD-1

## Files

This directory contains five files:
1. mkpsf_5ggr.tcl -- VMD/psfgen script that creates the first vacuum psf/pdb pair.  It uses a monte-carlo-based loop model-builder to build in the missing residues. 
2. my_5ggr_vac.namd -- NAMD configuration file used to relax the "guessed" coordinates resulting from step 1.
3. my_5ggr_solv.tcl -- VMD script that uses solvate and autoionize to generate a neutralized, solvated MD system using the coordinates from step 2 as input.
4. my_5ggr_colvars_op.inp -- colvars input file that defines collective variables that allow for center-of-mass restraint and an orientational restraint to keep the pseudo C2v axis of the Fab along z.
5. my_5ggr_solv.namd -- NAMD configuration file that performs a minimization and short MD of the raw solvated system; uses colvars module input file from 4.

## Instructions

Make sure PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local NAMD2 executable.  For me, these are /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

1. Download 5ggr.pdb to a clean directory

> wget www.rcsb.org/pdb/files/5ggr.pdb

2. Use VMD in text mode to generate the psf/pdb

> vmd -dispdev text -e $PSFGEN_BASEDIR/5ggr/mkpsf_5ggr.tcl

3. Run NAMD to relax bonds and guessed-in atoms

> ln -s $PSFGEN_BASEDIR/5ggr/my_5ggr_vac.namd .

> $CHARMRUN +p1 $NAMD2 my_5ggr_vac.namd > vac.log

4. Use VMD to solvate and neutralize the output of step 3

> vmd -dispdev text -e $PSFGEN_BASEDIR/5ggr/my_5ggr_solv.tcl

5. Run NAMD to minimize and shake out the solvated system

> ln -s $PSFGEN_BASEDIR/5ggr/my_5ggr_trimer_solv.namd .

> $CHARMRUN +p16 $NAMD2 my_5ggr_solv.namd > solv.log

2017, Cameron F Abrams
