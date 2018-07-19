# Some auxiliary scripts 

1. do_test.sh -- a Bash script that executes the general workflow to generate
   a solvated system.  Execute as

   > $PSFGEN_BASEDIR/scripts/do_test.sh -pdb 1abc

   where 1abc is replace with the name of a PDB file supported.  Do this in a clean directory.  Note that many of the systems in this repository have the own custom master scripts.

2. vmdrc.tcl -- this is my .vmdrc file (actually, my real .vmdrc just sources this file). Some helpful functions.

(c) 2017, Cameron F. Abrams
cfa22@drexel.edu
