# 3TGQ -- unliganded core monomeric HIV-1 gp120

[View PDB Entry](http://www.rcsb.org/structure/3TGQ)

[Read the source paper for this structure](http://dx.doi.org/10.1073/pnas.1112391109)

This entry uses the `cfapdbparser.py` python program to generate the Tcl `psfgen` script, with the whole build 
driven by `do_py.sh`.  The 3tgq entry specifies six unique biological assemblies, any one of which can be selected
as the build target (along with the asymmetric unit).  

Below are all the commands that are necessary.  Replace `<N>` with the number of the biological assembly you want.  You
may want to read the headers in the PDB file to decide; briefly, the asymmetric unit is a grouping of four distinct
gp120's, biological assemblies 1-4 are each individual gp120 alone (chains A, B, C, and D, respectively), while
assembles 5 and 6 are particular dimers.  Per usual, if no `-ba <N>` argument is passed to `cfapyparser.py`,
a system using the asymmetric unit is built.

```bash
> mkdir my_3tgq
> cd my_3tgq
> $PSFGEN_BASEDIR/scripts/do_py.sh -pyparserargs "-ba <N> -smdclose" -solv-stage-steps 10,200,100,1000,20000 -pdb 3tgq
```

There are a few short unresolved runs of residues that are modeled in by default, and the `-smdclose` directive
specifies that a step of steered MD is used to close the loops that result from this modeling in.  MD NPT equilibration 
of the solvated system here proceeds through five progressively longer stages to permit recalculation of the patch
grid.

2017-2021, Cameron F Abrams, cfa22@drexel.edu
