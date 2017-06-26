# 3TGQ -- unliganded core monomeric HIV-1 gp120

## Files

This directory contains four files:
1. mkpsf_3tgq.tcl -- this is the VMD/psfgen script that creates the first vacuum psf/pdb pair.  It uses a monte-carlo-based loop model-builder to build in the missing residues.
2. my_3tgq_vac.namd -- this is a NAMD configuration file used to relax the "guessed" coordinates resulting from mkpsf_3tgq.tcl.
3. solv.tcl -- this is a VMD script that uses solvate and autoionize to generate a neutralized, solvated MD system.
4. solv.namd -- this is a NAMD configuration file that performs a minimization and short MD of the raw solvated system.

## Instructions

If you have cloned this repository, then make sure PSFGEN_BASEDIR resolves to the root directory of your local copy.  If you did not
clone the repository, you will have to figure this part out on your own.

1. Download 3tgq.pdb to a clean directory

> wget http://www.rcsb.org/files/pdb/3tgq.pdb

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
