# 1HHP -- apo dimeric HIV-1 protease

## Files

This directory contains five files:
1. mkpsf_1hhp.tcl -- this is the VMD/psfgen script that creates the first vacuum psf/pdb pair.  It uses a monte-carlo-based loop model-builder to build in the missing residues.  It includes all glycans present in the PDB.
2. my_1hhp_vac.namd -- this is a NAMD configuration file used to relax the "guessed" coordinates resulting from step 1.
3. my_1hhp_solv.tcl -- this is a VMD script that uses solvate and autoionize to generate a neutralized, solvated MD system using the coordinates from step 2 as input.
4. my_1hhp_colvars_op.inp -- a colvars input file that defines collective variables that allow for center-of-mass restraint and an orientational restraint.
5. my_1hhp_solv.namd -- this is a NAMD configuration file that performs a minimization and short MD of the raw solvated system; uses colvars module input file from 4.

## Instructions

If you have cloned this repository, then make sure PSFGEN_BASEDIR resolves to the root directory of your local copy.  If you did not
clone the repository, you will have to figure this part out on your own.  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local NAMD2 executable.  For me, these are /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

1. Download 1hhp.pdb to a clean directory

> wget http://www.rcsb.org/pdb/files/1hhp.pdb

2. Use VMD in text mode to generate the psf/pdb

> vmd -dispdev text -e $PSFGEN_BASEDIR/1hhp/mkpsf_1hhp.tcl

3. Run NAMD to relax bonds and guessed-in atoms

> ln -s $PSFGEN_BASEDIR/1hhp/my_1hhp_vac.namd .

> $CHARMRUN +p1 $NAMD2 my_1hhp_vac.namd > vac.log

4. Use VMD to solvate and neutralize the output of step 3

> vmd -dispdev text -e $PSFGEN_BASEDIR/1hhp/my_1hhp_solv.tcl

5. Run NAMD to minimize and shake out the solvated system

> ln -s $PSFGEN_BASEDIR/1hhp/my_1hhp_trimer_solv.namd .

> $CHARMRUN +p16 $NAMD2 my_1hhp_solv.namd > solv.log

2017, Cameron F Abrams
