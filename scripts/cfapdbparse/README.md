# cfapdbparse.py

This python script automates the process of building an input script for `psfgen` using as much information as possible in a published PDB or mmCIF file.  

Cameron F Abrams -- cfa22@drexel.edu

2020-2021


### What it does

`cfapdbparse.py` parses most of the PDB/mmCIF records necessary to build a complete psfgen input script.  It uses SEQRES records to build chains, the REMARK 465 records to indicate missing residues, and the LINK and SSBOND records to include the proper bond patches.  It will generate replica protomers using BIOMT records to generate multimers.  It can generate the proper commands to mutate residues and cleave chains.

`cfapdbparse.py` is still a work in progress and has not been widely tested on random PDB entries.  Don't use it for any systems with ligands or co-factors (yet).  So far, the only HET's recognized are glycans.

A generic workflow might look like this:

```
> python cfadpbparse.py -pdb ####.pdb
> vmd -dispdev text mkpsf.tcl
> charmrun +pX namd2 my_####_vac.namd
```

The first command generates the `mkpsf.tcl` psfgen script; the second runs it to generate the initial psf/pdb pair, and the third runs the warmup simulation prior to solvation.
This workflow is expanded an automated in the `$PSFGEN_BASEDIR/scripts/do_py.sh` bash script.

By default, `cfapdbparse.py` assumes the input PDB file is an official entry in the PDB.  However, if it detects the keyword CHARMM in the KEYWDS record, it will then assume the PDB is in CHARMM format (in particular with regard to atom names and residue names).  This allows for chained processing, where a preliminary system is created and warmed-up, and then subsequently reprocessed through psfgen to add some other components or features.  This feature of `cfapdbparse.py` is functional, but the filename generation is not optimized.

### Specifying modifications

Modifications such as point mutations, point deletions, and grafting of glycans can be specified at the command line individually, or collectively in a `modsfile`.  The options `-mut` `-delete` and `-graft` can each be supplied with one or more arguments that specify these modifications, respectively.  However, for ease of reproducibility, it is advisable to put all modifications into a `modsfile`.  This is a text file that is divided into sections, each beginning with a bracketed keyword and ending with a blank line, with one record per line.  For instance, a modfile that implements mutations and deletions of the SARS-CoV-2 spike based on the 6vxx PDB entry in order to implement the Delta variant sequence might look like this:

```
[title]
SARS-CoV-2 Delta S

[description]
Modifications to 6vxx PDB structure to generate the B.1.617.2 Delta variant, fully glycosylated
Step 1: Pre-cleavage modifications

[mutations]
A_L452R
A_T478K
A_D614G
A_P681R
A_D950N
B_L452R
B_T478K
B_D614G
B_P681R
B_D950N
C_L452R
C_T478K
C_D614G
C_P681R
C_D950N

[deletions]
A_F157
A_R158
B_F157
B_R158
C_F157
C_R158
```

### Requirements

Works with python 3.x, with `numpy` and `pycifrw` packages required.
