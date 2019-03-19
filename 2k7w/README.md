# 2K7W -- BAX with bound BIM SAHB peptide

## Instructions

Make sure PSFGEN_BASEDIR resolves to the root directory of your local copy of this repository (mine is ${HOME}/research/psfgen).  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local namd2 executable.  For me, these are ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and ${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

Use the generate bash script `do_test.sh` to generate a run ready system in a clean directory:
```
mkdir my_bax
cd my_bax
${PSFGEN_BASEDIR}/scripts/do_test.sh -pdb 2k7w [-psfgen_args -frm # [-nobim] [-P168G] [-cisP168] ]
```

flag | default | description
---|---|---
`-frm` | (last frame) | signify which frame from the NMR dataset to use
`-nobim` | n/a | exclude BIM SAHB peptide
`-P168G` | n/a | mutate proline 168 to glycine
`-cisP168` | n/a | put the 167-168 peptide bond omega angle into cis

Note that the `-cisP168` flag triggers a second stage of vacuum MD prior to solvation.  In this stage, the colvars module is used to force the 167-168 peptide bond into cis.

2019, Cameron F Abrams
