# 4H8W -- core monomeric HIV-1 clade A/E gp120

## Files

This directory contains five files:
1. mkpsf_4h8w.tcl --  VMD/psfgen script that creates the first vacuum psf/pdb pair.  It uses a monte-carlo-based loop model-builder to build in the missing residues.
2. my_4h8w_vac.namd --  NAMD configuration file used to relax the "guessed" coordinates resulting from step 1, which iclude all hydrogens and all atoms on model-built loops.
3. my_4h8w_solv.tcl -- VMD script that uses solvate and autoionize to generate a neutralized, solvated MD system.
4. my_4h8w_colvars_op.inp -- a colvars input file that defines collective variables that allow for center-of-mass restraint and orientational restraint.
5. my_4h8w_solv.namd -- a NAMD configuration file that performs a minimization and short MD of the raw solvated system.

## Quick Instructions

The file do_test.sh in the scripts/directory is a Bash script that will generate a run-ready, relaxed, solvated system if issued in a clean directory.  The steps it performs are detailed in the following section.  In a clean directory, issue:

> $PSFGEN_BASEDIR/scripts/do_test.sh -pdb 4h8w

NOTE:  To revert the serine at 375 back to a histidine that is considered WT for clade A/E HIV-1, issue the command as

> $PSFGEN_BASEDIR/scripts/do_test.sh -pdb 4h8w -psfgen_args S375H

(Other mutations will soon be supported -- this is just a quick one)

## Detailed Instructions

Make sure PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local NAMD2 executable.  For me, these are ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

1. Download 4h8w.pdb to a clean directory

> wget http://www.rcsb.org/pdb/files/4h8w.pdb

2. Use VMD in text mode to generate the psf/pdb

> vmd -dispdev text -e $PSFGEN_BASEDIR/4h8w/mkpsf_4h8w.tcl

3. Run NAMD to relax bonds and guessed-in atoms

> ln -s $PSFGEN_BASEDIR/4h8w/my_4h8w_vac.namd .

> $CHARMRUN +p1 $NAMD2 my_4h8w_vac.namd > vac.log

4. Use VMD to solvate and neutralize the output of step 3

> vmd -dispdev text -e $PSFGEN_BASEDIR/4h8w/solv.tcl

5. Run NAMD to minimize and shake out the solvated system

> ln -s $PSFGEN_BASEDIR/4h8w/solv.namd .

> $CHARMRUN +p16 $NAMD2 solv.namd > solv.log


2017, Cameron F Abrams
