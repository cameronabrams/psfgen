# 7JJI -- Soluble, stabilized SARS-Cov-2 S trimeric spike, closed, with linoleic acid and polysorbate-80 ligands, C3 symmetric

[View PDB Entry](http://www.rcsb.org/structure/7JJI)

[Read the source paper for this structure](http://dx.doi.org/10.1101/2020.08.06.234674)

## Files provided here

`grafts.inp`:  A list of glycans to be grafted onto the main molecule.  Each item in the list is of the form
```
PDB.pdb,X:#-#,Y:#,#
```
where `PDB.pdb` is the name of pdb file in the RCSB that contains the glycan, and `X:#-#` indicates the chain and resid span for that glycan in that file.  `Y:#` designates the target chain and resid onto which the graft is placed.  Note that the target residue and the **first** residue in the graft must have heavy-atom order/type congruency (should be same resname).  Finally, the last number is a desired resid offset; this number is added to the **target** resid and all subsequent residues in the graft are "downstream" by resid of that number.  With multiple grafts, it is usually a good idea to set the offsets so they are at least 100 apart and 100 greater then the highest resid in the existing target molecule.

`crot.inp`: A list of explicit bond rotations that move glycan fragments out of steric overlap with other fragments.

## Instructions

This workflow generates a solvated, cleaved, fully glycosylated SARS-CoV-2 S spike ectodomain trimer based on the 7jji PDB entry.  It uses the `cfapdbparser.py` package and the general driver `do_py.sh`.   Make sure your environment variable PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local NAMD2 executable.  (For me, these are /home/cfa/namd/NAMD_2.14_Source/Linux-x86_64-g++/charmrun and /home/cfa/namd/NAMD_2.14_Source/Linux-x86_64-g++/namd2.

```
$ mkdir 7jji
$ cd 7jji
$ $PSFGEN_BASEDIR/scripts/do_py.sh -pyparser-args "-grafile $PSFGEN_BASEDIR/7jji/grafts.inp -crotfile $PSFGEN_BASEDIR/7jji/crot.inp -loctopo vcg-paramchem-ic.str -locparam vcg-paramchem-ic.str -charmmtopo stream/lipid/toppar_all36_lipid_detergent.str -charmmparam stream/lipid/toppar_all36_lipid_detergent.str -rem -smdheal" -pyparser-args "-clv A685 -clv B685 -clv C685 -loctopo vcg-paramchem-ic.str -locparam vcg-paramchem-ic.str -charmmtopo stream/lipid/toppar_all36_lipid_detergent.str -charmmparam stream/lipid/toppar_all36_lipid_detergent.str" -solv-stage-steps 100,200,400,800,20000 -temperature 310 -pdb 7jji -pdb 2wah -pdb 4byh -pdb 4b7i
```

The `do_py.sh` script executes a series of tasks, beginning with downloading the required PDB file from the RCSB (if needed), then passing through a sequence of parse/psfgen/relax cycles to generate a complete vacuum structure, followed by solvation via psfgen, and finally through as series of solvated relaxations via NPT MD.  

In this particular case, the driver runs two parser instances in series, indicated by the two `-pyparser-args` arguments.  The first will add missing loops and graft on glycans, and the second executes the cleavages at the furin sites.  Also in the first round, the `-rem` switch reverses all ENGINEERED MUTATIONS listed in the source PDB file.  The `-smdheal` switch enables the steered-MD loop relaxation step (preferred over loop Monte Carlo).  The additional topology and parameter files take care of the linoleic acid (in the standard CHARMM detergent params) and the `vcg` files take care of polysorbate-80, which I parameterized using paramchem.  Five stages of solvated equilibration are requested (which helps with patch-grid errors as the box size equilibrates).  The 2wah, 4byh, and 4b7i PDB entries contain large glycans that are grafted according to the graft records in `grafts.inp`.

To make an uncleaved S trimer, simply omit the second `-pyparser-args` switch.

The glycans are assigned according to [Watanabe et al.](https://science.sciencemag.org/content/369/6501/330) with glycans classified as "complex" represented by the glycan in 4byh, "hybrid" with 4b7i, and oligomannose with 2wah.

2017-2020, Cameron F Abrams  cfa22@drexel.edu
