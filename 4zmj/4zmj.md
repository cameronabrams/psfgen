# 4ZMJ -- unliganded SOSIP trimer with glycans

## Files

This directory contains five files:
1. mkpsf_4zmj_trimer.tcl -- this is the VMD/psfgen script that creates the first vacuum psf/pdb pair.  It uses a monte-carlo-based loop model-builder to build in the missing residues.  It includes all glycans present in the PDB entry and performs symmetry replication as instructed to generate the full trimer.
2. my_4zmj_trimer_vac.namd -- this is a NAMD configuration file used to relax the "guessed" coordinates resulting from mkpsf_4zmj_trimer.tcl.
3. my_4zmj_trimer_solv.tcl -- this is a VMD script that uses solvate and autoionize to generate a neutralized, solvated MD system using the coordinates from step 2 as input.
4. my_4zmj_trimer_colvars_op.inp -- a colvars input file that defines collective variables that allow for center-of-mass restraint and an orientational restraint to keep the C3v axis along z.
5. my_4zmj_trimer_solv.namd -- this is a NAMD configuration file that performs a minimization and short MD of the raw solvated system; uses colvars module input file from 4.

## Instructions

If you have cloned this repository, then make sure PSFGEN_BASEDIR resolves to the root directory of your local copy.  If you did not
clone the repository, you will have to figure this part out on your own.  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local NAMD2 executable.  For me, these are /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

1. Download 4zmj.pdb to a clean directory

> wget http://www.rcsb.org/files/pdb/4zmj.pdb

2. Use VMD in text mode to generate the psf/pdb

> vmd -dispdev text -e $PSFGEN_BASEDIR/4zmj/mkpsf_4zmj_trimer.tcl

3. Run NAMD to relax bonds and guessed-in atoms

> ln -s $PSFGEN_BASEDIR/4zmj/my_4zmj_trimer_vac.namd .

> $CHARMRUN +p1 $NAMD2 my_4zmj_trimer_vac.namd > vac.log

4. Use VMD to solvate and neutralize the output of step 3

> vmd -dispdev text -e $PSFGEN_BASEDIR/4zmj/my_4zmj_trimer_solv.tcl

5. Run NAMD to minimize and shake out the solvated system

> ln -s $PSFGEN_BASEDIR/4zmj/my_4zmj_trimer_solv.namd .

> $CHARMRUN +p16 $NAMD2 my_4zmj_trimer_solv.namd > solv.log

2017, Cameron F Abrams
