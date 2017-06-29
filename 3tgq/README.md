# 3TGQ -- unliganded core monomeric HIV-1 gp120

## Files

This directory contains five files:
1. mkpsf_3tgq.tcl --  VMD/psfgen script that creates the first vacuum psf/pdb pair.  It uses a monte-carlo-based loop model-builder to build in the missing residues.
2. my_3tgq_vac.namd --  NAMD configuration file used to relax the "guessed" coordinates resulting from step 1, which iclude all hydrogens and all atoms on model-built loops.
3. my_3tgq_solv.tcl -- VMD script that uses solvate and autoionize to generate a neutralized, solvated MD system.
4. my_3tgq_colvars_op.inp -- a colvars input file that defines collective variables that allow for center-of-mass restraint and orientational restraint.
5. my_3tgq_solv.namd -- a NAMD configuration file that performs a minimization and short MD of the raw solvated system.

The file do_test.sh in the scripts/directory is a Bash script that performs the sequence of commands detailed in the instructions below.

## Instructions

Make sure PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local NAMD2 executable.  For me, these are ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

1. Download 3tgq.pdb to a clean directory

> wget http://www.rcsb.org/pdb/files/3tgq.pdb

2. Use VMD in text mode to generate the psf/pdb

> vmd -dispdev text -e $PSFGEN_BASEDIR/3tgq/mkpsf_3tgq.tcl

3. Run NAMD to relax bonds and guessed-in atoms

> ln -s $PSFGEN_BASEDIR/3tgq/my_3tgq_vac.namd .

> $CHARMRUN +p1 $NAMD2 my_3tgq_vac.namd > vac.log

4. Use VMD to solvate and neutralize the output of step 3

> vmd -dispdev text -e $PSFGEN_BASEDIR/3tgq/solv.tcl

5. Run NAMD to minimize and shake out the solvated system

> ln -s $PSFGEN_BASEDIR/3tgq/solv.namd .

> $CHARMRUN +p16 $NAMD2 solv.namd > solv.log


2017, Cameron F Abrams
