# 5F4P -- core monomeric HIV-1 gp120 with smCD4mc BNM-III-170 bound

## About this structure

Small-molecule CD4-mimetic compounds [1] bind to HIV-1 gp120 in the 
Phe43 cavity and induce a conformational response in HIV-1 env
that mimics that induced by CD4 binding.  In this structure, the
Hendrickson group showed that, like other smCD4mc's, the particular
one BNM-III-170 binds by placing its chlorofluorophenyl deep inside
the cavity [2].  Here, I have used CGenFF to determine atom types
for the BNM molecule, which are stored in charmm/bnm.str.  The only
manual change was to increase the kappa value on the C=O-C=O torsion
on the oxalamide to 4 to keep the oxalamide planar.

## Instructions (bash)

```
$ export PSFGEN_BASEDIR=/path/to/psfgen/repository/on/your/system
$ mkdir my_5f4p
$ cd my_5f4p
$ $PSFGEN_BASEDIR/5f4p/do_test.sh
```

## References

1. "Structure-based design, synthesis and validation of CD4-mimetic small molecule inhibitors of HIV-1 entry: conversion of a viral entry agonist to an antagonist.", Courter JR, Madani N, Sodroski J, Schoen A, Freire E, Kwong PD, Hendrickson WA, Chaiken IM, LaLonde JM, Smith AB III.  _Acc Chem Res._ 2014;47(4):1228-37. doi: 10.1021/ar4002735.

2. "Small-Molecule CD4-Mimics: Structure-Based Optimization of HIV-1 Entry Inhibition.", Melillo, B., Liang, S., Park, J., Schoen, A., Courter, J.R., LaLonde, J.M., Wendler, D.J., Princiotto, A.M., Seaman, M.S., Freire, E., Sodroski, J., Madani, N., Hendrickson, W.A., Smith, A.B. III. _Acs Med.Chem.Lett._ 2016;7:330-334

2018, Cameron F Abrams
