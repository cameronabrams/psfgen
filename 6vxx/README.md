# 6VXX -- Soluble, stabilized SARS-Cov-2 S trimeric spike, closed

[View PDB Entry](http://www.rcsb.org/structure/6VXX)

[Read the source paper](http://dx.doi.org/10.1016/j.cell.2020.02.058)

## Instructions

This workflow for generating a solvated, cleaved, fully glycosylated spike ectodomain trimer uses the new `cfapdbparser.py` package and the general driver `do_py.sh`.   A working python distribution is assumed.  Make sure your environment variable PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local NAMD2 executable.  (For me, these are /home/cfa/namd/NAMD_2.13_Source/Linux-x86_64-g++/charmrun and /home/cfa/namd/NAMD_2.13_Source/Linux-x86_64-g++/namd2.

```
$ mkdir 6vxx
$ cd 6vxx
$ $PSFGEN_BASEDIR/scripts/do_py.sh -pyparser-args '-graftfile $PSFGEN_BASEDIR/6vxx/grafts.in -rmi' -pyparser-args '-clv A685S -clv B685T -clv C685U' -solv-stage-steps 100,200,400,800,20000 -temperature 310 -pdb 6vxx -pdb 2wah -pdb 4byh -pdb 4b7i
```

This tells the driver to run two parser instances in series.  The first will add missing loops and graft on glycans, and the second executes the cleavage at the furin site.  Five stages of solvate equilibration are requested (which helps with patch-grid errors as the box size equilibrates.  The 2wah, 4yh, and 4b7i PDB entries contain large glycans that are grafted according to the graft records in `grafts.in`.

The glycans are assigned according to [Watanabe et al.](https://science.sciencemag.org/content/369/6501/330).

2017-2020, Cameron F Abrams
