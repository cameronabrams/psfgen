# 4ZMJ -- Soluble, stabilized HIV-1 Env trimeric spike ectodomain, closed

[View PDB Entry](http://www.rcsb.org/structure/4ZMJ)

[Read the source paper for this structure](http://dx.doi.org/10.1038/nsmb.3051)

## Files provided here

`crot.inp`:  A list of phi, psi rotations used to help insert model-built loops

## Instructions

This workflow generates a solvated, cleaved, partially glycosylated HIV-1 Env ectodomain trimer based on the 4zmj PDB entry.  It uses the `cfapdbparser.py` package and the general driver `do_py.sh`.   Make sure your environment variable PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local NAMD2 executable.  (For me, these are /home/cfa/namd/NAMD_2.13_Source/Linux-x86_64-g++/charmrun and /home/cfa/namd/NAMD_2.13_Source/Linux-x86_64-g++/namd2.

```
$ mkdir my_4zmj
$ cd my_4zmj
$ $PSFGEN_BASEDIR/scripts/do_py.sh -pyparser-args '-crotfile $PSFGEN_BASEDIR/4zmj/crot.inp -rmi' -solv-stage-steps 100,200,400,800,20000 -temperature 310 -pdb 4zmj
```

The `do_py.sh` script executes a series of tasks, beginning with downloading the required PDB file from the RCSB (if needed), then passing through a sequence of parse/psfgen/relax cycles to generate a complete vacuum structure, followed by solvation via psfgen, and finally through as series of solvated relaxations via NPT MD.  

In this particular case, the driver runs one parser instance, indicated by the single `-pyparser-args` argument.  The will add missing loops and with some minor torsion adjustments to avoid clashes. Five stages of solvated equilibration are requested (which helps with patch-grid errors as the box size equilibrates.  

2017-2020, Cameron F Abrams
