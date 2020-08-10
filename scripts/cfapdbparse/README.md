# cfapdbparse.py

This python script automates the process of building an input script for psfgen using as much information as possible in the published PDB file.  

Cameron F Abrams -- cfa22@drexel.edu

2020


### How it works:

`cfapdbparse` parses most of the PDB records necessary to build a complete psfgen input script.  It uses SEQRES records to build chains, the REMARK 465 records to indicate missing residues, and the LINK and SSBOND records to include the proper bond patches.  It can generate the proper commands to mutate residues and cleave chains.

`cfapdbparse` is still a work in progress and has not been widely tested on random PDB entries.  Don't use it for any systems with ligands or co-factors (yet).  So far, the only HET's recognized are glycans.

A generic workflow might look like this:

```
> python cfadpbparse.py -pdb ####.pdb
> vmd -dispdev text mkpsf.tcl
> charmrun +pX namd2 my_####_vac.namd
```

The first command generates the `mkpsf.tcl` psfgen script; the second runs it to generate the initial psf/pdb pair, and the third runs the warmup simulation prior to solvation.

By default, `cfapdbparse` assumes the input PDB file is an official entry in the PDB.  However, if it detects the keyword CHARMM in the KEYWDS record, it will then assume the PDB is in CHARMM format (in particular with regard to atom names and residue names).  This allows for chained processing, where a preliminary system is created and warmed-up, and then subsequently reprocessed through psfgen to add some other components or features.  This feature of `cfapdbparse` is functional, but the filename generation is not optimized.
