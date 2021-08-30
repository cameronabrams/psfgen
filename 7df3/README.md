# 7DF3 -- Soluble, stabilized SARS-Cov-2 S trimeric spike, closed

[View PDB Entry](http://www.rcsb.org/structure/7DF3)

[Read the source paper for this structure](http://dx.doi.org/10.1126/sciadv.abe5575)

## Files provided here

1. `glycans.mod`:  A specially-formatted "modsfile" lising glycans to be grafted onto the main molecule; generated using `glycan-extract.py` (included here for completeness).
2. `delta.mod`: A modsfile that specifies mutations and deletions to generate a B.1.617.2 Delta variant.

Each item in a `[grafts]` stanza of a modsfile is of the form
    ```
    PDB.pdb,X:#-#,Y:#,#
    ```
where `PDB.pdb` is the name of pdb file in the RCSB that contains the glycan, and `X:#-#` indicates the chain and resid span for that glycan in that file.  `Y:#` designates the target chain and resid onto which the graft is placed, typically a stem NAG linked to an asparagine's ND2 atom.  Note that the target residue and the **first** residue in the graft must have heavy-atom order/type congruency (should be same resname).  Finally, the last number is a desired resid offset; this number is added to the **target** resid and all subsequent residues in the graft are "downstream" by resid of that number.  With multiple grafts, it is usually a good idea to set the offsets so they are at least 100 apart and 100 greater then the highest resid in the existing target molecule.

## Instructions

This workflow generates a solvated, cleaved, fully glycosylated SARS-CoV-2 S spike ectodomain trimer based on the 6vxx PDB entry.  It uses the `cfapdbparser.py` package and the general driver `do_py.sh`.   Make sure your environment variable PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local NAMD2 executable.  (For me, these are /home/cfa/namd/NAMD_2.14_Source/Linux-x86_64-g++/charmrun and /home/cfa/namd/NAMD_2.14_Source/Linux-x86_64-g++/namd2.

```bash
> cp -r $PSFGEN_BASEDIR/7df3 .
> cd 7df3
> $PSFGEN_BASEDIR/scripts/do_py.sh -pyparser-args "-modsfile glycans.mod -smdclose" \
                                   -pyparser-args "-clv A685 B685 C685" \
                                   -solv-stage-steps 100,200,400,800,20000 \
                                   -temperature 310 -pdb 7df3 2wah 4byh 4b7i
```

The `do_py.sh` script executes a series of tasks, beginning with downloading the required PDB file from the RCSB (if needed), then passing through a sequence of parse/psfgen/relax cycles to generate a complete vacuum structure, followed by solvation via psfgen, and finally through as series of solvated relaxations via NPT MD.  

In this particular case, the driver runs two parser instances in series, indicated by the two `-pyparser-args` arguments.  The first will add missing loops and graft on glycans, and the second executes the cleavages at the furin sites.  The `-smdclose` switch indicates that missing loops that are added by psfgen will be "closed" using steered MD.  Five stages of solvated equilibration are requested (which helps with patch-grid errors as the box size equilibrates).  The 2wah, 4yh, and 4b7i PDB entries contain large glycans that are grafted according to the lines in the `[grafts]` stanza in the modsfile `grafts.inp`.

To make an uncleaved S trimer, simply omit the second `-pyparser-args` switch.

The glycans are assigned according to [Watanabe et al.](https://science.sciencemag.org/content/369/6501/330) with glycans classified as "complex" represented by the glycan in 4byh, "hybrid" with 4b7i, and oligomannose with 2wah.

Any desired point mutations can be specified in the `[mutations]` section of the modsfile, and any point deletions can be specified in the `[deletions]` section.  
The format of a point mutation is `C_OrrrN` where `C` is a chain ID, `O` is the original one-letter residue name, `rrr` is the residue sequence number, and `N` is the desired mutant one-letter residue name.  Point-deletions are specified like mutations except with no `N`.  Each line in any modsfile section can contain only one specification.  

To build a spike that conforms to the Delta variant sequence:
```bash
> cp -r $PSFGEN_BASEDIR/7df3 .
> cd 7df3
> $PSFGEN_BASEDIR/scripts/do_py.sh -pyparser-args "-modsfile glycans.mod delta.mod -smdclose" \
                                   -pyparser-args "-clv A685 B685 C685" \
                                   -solv-stage-steps 100,200,400,800,20000 \
                                   -temperature 310 -pdb 7df3 2wah 4byh 4b7i
```

2017-2021, Cameron F Abrams  cfa22@drexel.edu
