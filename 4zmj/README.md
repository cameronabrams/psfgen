# 4ZMJ -- unliganded SOSIP trimer with glycans; option to model-in MPER's

## Files

This directory contains five files:
1. mkpsf_4zmj.tcl -- this is the VMD/psfgen script that creates the first vacuum psf/pdb pair.  It uses a monte-carlo-based loop model-builder to build in the missing residues.  It includes all glycans present in the PDB entry and performs symmetry replication as instructed to generate the full trimer.
2. my_4zmj_vac.namd -- this is a NAMD configuration file used to relax the "guessed" coordinates resulting from mkpsf_4zmj.tcl.
3. my_4zmj_solv.tcl -- this is a VMD script that uses solvate and autoionize to generate a neutralized, solvated MD system using the coordinates from step 2 as input.
4. my_4zmj_colvars_op.inp -- a colvars input file that defines collective variables that allow for center-of-mass restraint and an orientational restraint to keep the C3v axis along z.
5. my_4zmj_solv.namd -- this is a NAMD configuration file that performs a minimization and short MD of the raw solvated system; uses colvars module input file from 4.

## Instructions

```
mkdir my_4zmj
cd my_4zmj
$PSFGEN_BASEDIR/scripts/do_test.sh -pdb 4zmj [-psfgen_args -seed # -mper-extend]
```

The `-seed` flag allows the user to set the random-number generator seed.  The flag `-mper-extend` instructs the script to grow in the MPER residues on each gp41 as an alpha-helix, out to residue 682.

2017-2018, Cameron F Abrams
cfa22@drexel.edu
