# 5X58 -- Soluble, stabilized SARS-Cov-1 S trimeric spike, closed (conformation 1)

[View PDB Entry](http://www.rcsb.org/structure/5X5B)

[Read the source paper for this structure](https://www.nature.com/articles/ncomms15092)

## Files provided here

  * `userlink.inp`:  A list PDB-format LINK records that seem to be missing from the RCSB PDB.

  * `glycangrafts.inp`: A list of glycan graft records, each of the form
```
PDB.pdb,X:#-#,Y:#,#
```
where `PDB.pdb` is the name of pdb file in the RCSB that contains the glycan, and `X:#-#` indicates the chain and resid span for that glycan in that file.  `Y:#` designates the target chain and resid onto which the graft is placed.  Note that the target residue and the **first** residue in the graft must have heavy-atom order/type congruency (should be same resname).  Finally, the last number is a desired resid offset; this number is added to the **target** resid and all subsequent residues in the graft are "downstream" by resid of that number.  With multiple grafts, it is usually a good idea to set the offsets so they are at least 100 apart and 100 greater then the highest resid in the existing target molecule.

## Instructions

This workflow generates a solvated, fully glycosylated SARS-CoV-1 S spike ectodomain trimer based on the 5x58 PDB entry.  It uses the `cfapdbparser.py` package and the general driver `do_py.sh`.   Make sure your environment variable PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local NAMD2 executable.  (For me, these are /home/cfa/namd/NAMD_2.14_Source/Linux-x86_64-g++/charmrun and /home/cfa/namd/NAMD_2.14_Source/Linux-x86_64-g++/namd2.

```
$ mkdir 5x58
$ cd 5x58
$ $PSFGEN_BASEDIR/scripts/do_py.sh -pyparser-args "-grafile $PSFGEN_BASEDIR/5x58/glycangrafts.inp -linkfile $PSFGEN_BASEDIR/5x58/userlink.inp -rmi" -solv-stage-steps 100,200,400,800,20000 -temperature 310 -pdb 5x58 -pdb 2wah -pdb 4byh -pdb 4b7i
```

The `do_py.sh` script executes a series of tasks, beginning with downloading the required PDB file from the RCSB (if needed), then passing through a sequence of parse/psfgen/relax cycles to generate a complete vacuum structure, followed by solvation via psfgen, and finally through as series of solvated relaxations via NPT MD.  

In this particular case, the driver runs a single parser instance, indicated by the `-pyparser-args` argument.  This will add missing loops and graft on glycans. Five stages of solvated equilibration are requested (which helps with patch-grid errors as the box size equilibrates).  The 2wah, 4yh, and 4b7i PDB entries contain large glycans that are grafted according to the graft records in `grafts.inp`.

The glycans are assigned according to [Walls et al.](https://doi.org/10.1016/j.cell.2018.12.028) with glycans classified as "complex" represented by the glycan in 4byh, "hybrid" with 4b7i, and oligomannose with 2wah.  Of particular note: Walls et al. only detected a NAG at N158, and a NAG is present at 158 in 5x58, so we do not graft onto that.  Also, N318 and N602 only showed glycans in their EM structures, not in the LC-MS/MS data; we assume these are complex.  We graft the glycan they detect at N1053 to the 5x58 N1056.  We currently do not graft any glycan to N589, or any N's above and including 1080.

2017-2020, Cameron F Abrams  cfa22@drexel.edu
