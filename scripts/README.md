# Some auxiliary scripts 

1. `do_test.sh` -- a Bash script that executes the general workflow to generate
   a solvated system.  Execute as
   ```bash
   > $PSFGEN_BASEDIR/scripts/do_test.sh -pdb 1abc
   ```
   where 1abc is replace with the name of a PDB file supported.  Do this in a clean directory.  Note that many of the systems in this repository have the own custom master scripts.

2. `vmdrc.tcl` -- this is my `.vmdrc` file (actually, my real .vmdrc just sources this file). Some helpful functions.

3. `ringp.tcl` -- defines a function that is used to check for pierced rings in model-building

4. `tmd_prep.tcl` -- used for setting up TMD simulations

5. `solv.tcl` -- just solvates

6. `cfapdbparse` -- an experimental python program for automated system building; `do_py.sh` drives it

7. `tg.sh` and `tg.tcl` -- topogromacs driver for building gromacs topologies.  Note:  you **must** comment out the line `display update ui` in `topotools1.8/topogromacs.tcl` in the VMD plugins directory.

(c) 2017-2021, Cameron F. Abrams
cfa22@drexel.edu
