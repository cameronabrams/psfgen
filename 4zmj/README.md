# 4ZMJ -- Soluble, stabilized HIV-1 Env trimeric spike ectodomain, closed

[View PDB Entry](http://www.rcsb.org/structure/4ZMJ)

[Read the source paper for this structure](http://dx.doi.org/10.1038/nsmb.3051)

## Files provided here

`crot.inp`:  A list of phi, psi rotations used to help insert model-built loops

## Instructions

This workflow generates a solvated, cleaved, partially glycosylated HIV-1 Env ectodomain trimer based on the 4zmj PDB entry.  It uses the `cfapdbparser.py` package and the general driver `do_py.sh`.   Make sure your environment variable `PSFGEN_BASEDIR` resolves to the root directory of your local copy of this repository (mine is `${HOME}/research/psfgen`).  It is also assumed below that `CHARMRUN` resolves to your local `charmrun` executable and `NAMD2` resolves to your local `namd2` executable.

```
$ mkdir my_4zmj
$ cd my_4zmj
$  $PSFGEN_BASEDIR/scripts/do_py.sh -pyparser-args "-crotfile $PSFGEN_BASEDIR/4zmj/crot.inp -smdclose -ba 1 -fixconflicts -fixengineeredmutations" -solv-stage-steps 100,200,400,800,20000 -temperature 310 -pdb 4zmj -make-gromacs my_4zmj.top my_4zmj.pdb
```

The `do_py.sh` script executes a series of tasks, beginning with downloading the required PDB file from the RCSB (if needed), then passing through a sequence of parse/psfgen/relax tasks to generate a complete vacuum structure, followed by solvation via `psfgen`, and finally through as series of solvated relaxations via `NAMD`.  The `-smdclose` switch enables a steered MD task that connects the C-termini of modeled-in interior loops of residues indicated as "MISSING" with their correctly resolved C-terminal neighbor along each chain.  

The `-ba 1` switch selects Biological Assembly 1, which is the trimer.  Biological Assembly 0 is (always) the asymmetric unit, which in this case is a single protomer.  Leaving the `-ba` switch off forces a build of just the asymmetric unit, which is probably not what you want to do here.

In this particular case, the driver runs one parser instance, indicated by the single `-pyparser-args` argument.  The will add missing loops and with some minor torsion adjustments to avoid clashes (listed in the `crot.inp` file). Five stages of solvated equilibration are requested (which helps with patch-grid errors as the box size equilibrates.  

The `-fixconflicts` and `-fixengineeredmutations` take care of mutating residues explicitly listed in the PDB back to their 
original identities based on information in SEQADV records.  Conflicts refer to residue identities that do not match the
sequence identity for the construct; 4zmj has two conflicts, one of which leads to an N-linked glycosyation!  Engineered mutations are sequence changes introduced purposefully to aid in structure resolution; here, this refers to the SOSIP modifications, including the isoleucine-to-proline and the engineered disulfide linking gp120 to gp41.  With the directives above, all of 
these are mutated back to their wild-type identities.  It is important to note that the "mutations" these fixes introduce result in deletion of a glycan (chain D in the asymmetric unit) and deletion of a disulfide.

To introduce point mutations, use a -mut switch in the pyparser-args argument; e.g., for K117A in chain G:

```
$ $PSFGEN_BASEDIR/scripts/do_py.sh -pyparser-args "-crotfile $PSFGEN_BASEDIR/4zmj/crot.inp -smdclose -mut G_K117A -ba 1 -fixconflicts -fixengineeredmutations" -solv-stage-steps 100,200,400,800,20000 -temperature 310 -pdb 4zmj -make-gromacs my_4zmj.top my_4zmj.pdb
```
Because of the `BIOMT` transformations in the first biological assembly, this mutation will also be carried out on the symmetry-related chains of chain G.

2017-2021, Cameron F Abrams
