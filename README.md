# psfgen -- Advanced scripts for generating configurations and CHARMM36 topologies for use with NAMD

This repository contains psfgen scripts, TcL scripts for use in VMD, and some associated CHARMM36 topology and parameter files and NAMD config files used to generate _initial conditions_ for production MD simulations of various systems.  It should be helpful to anyone already familiar with using [psfgen](https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen) to build systems using CHARMM topologies who also needs access to some more advanced system-building capabilities than are available in the [psfgen tutorial](https://www.ks.uiuc.edu/Research/namd/tutorial/NCSA2002/hands-on).  It is strongly recommended that users have familiarity with the psfgen plugin via the excellent [user guide](https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/ug.pdf). Some of the features this repository provides beyond what psfgen can easily do are the following:

* support for rudimentary loop model-building to include residues missing in a PDB file but present in the crystallized protein sequence;
* support for glycans and non-covalently linked sugars and other ligands;
* integration between solvation and initial MD simulation config file (easy transfer of box size);
* support for down-puckered prolines;
* de-novo membrane-building using packmol

The _src_ directory contains files that can be sourced by psfgen scripts.  The _charmm_ directory contains some custom topologies and parameters derived from the July, 2016 charmm36 parameter set.  Other directory names indicate the PDB entry for which the files contained therein are applicable.

The repository is being updated continuously.  Issue `git pull` in your local copy often to keep it up-to-date.  The `master` branch should be functional.

## Requirements

1. NAMD v. 2.13.  Set environment variables CHARMRUN and NAMD2 to point to your system's `charmrun` and `namd2`; examples (in `~/.bashrc`):
```
export CHARMRUN=${HOME}/namd/NAMD_2.13_Source/Linux-x86_64-g++/charmrun
export NAMD2=${HOME}/namd/NAMD_2.13_Source/Linux-x86_64-g++/namd2
```
2. VMD v. 1.8.3
3. CHARMM36 topologies and parameters (toppar_c36_jul16.tgz is the version used here) unpacked in ${HOME}/charmm/toppar
4. packmol
5. tcl, tcl-devel, and swig
6. This repository cloned into a local directory and pointed to by the environment variable PSFGEN_BASEDIR; for example, in `~/.bashrc`,
```
export PSFGEN_BASEDIR=/home/myusername/psfgen
```
7. A locally compiled `bondstruct.so` module for the loop Monte Carlo procedures.  To build this:

```
$ cd $PSFGEN_BASEDIR
$ mkdir lib
$ cd src
$ make bondstruct.so
```

8. Add the following to your `~/.vmdrc` file:
```
source $PSFGEN_BASEDIR/scripts/vmdrc.tcl
```

## Structures

1. [3TGQ](http://www.rcsb.org/pdb/explore/explore.do?structureId=3tgq) -- unliganded core monomeric HIV-1 gp120, with NAG's;

2. [4ZMJ](http://www.rcsb.org/pdb/explore/explore.do?structureId=4zmj) -- unliganded soluble trimeric HIV-1 Env gp140, with glycans;

3. [5FUU](http://www.rcsb.org/pdb/explore/explore.do?structureId=5fuu) -- soluble, cleaved JR-FL HIV-1 trimer, with glycans; option to build Man9's on several sites

4. [1HHP](http://www.rcsb.org/pdb/explore/explore.do?structureID=1hhp) -- apo, dimeric HIV-1 protease;

5. [2MB5](http://www.rcsb.org/pdb/explore/explore.do?structureID=2mb5) -- carbonmonoxymyoglobin;

6. [1F7A](http://www.rcsb.org/pdb/explore/explore.do?structureID=1f7a) -- HIV-1 protease with bound substrate peptide;

7. [1HIW](http://www.rcsb.org/pdb/explore/explore.do?structureID=1hiw) -- HIV-1 matrix (MA) trimer;

8. [4H8W](http://www.rcsb.org/pdb/explore/explore.do?structureID=4h8w) -- HIV-1 gp120 core clade A/E, with options to implement any number of the following mutations: (1) S375H (to bring back to clade A/E WT), (2) H61Y, (3) Q105H, (4) V108I, or (5) NIK474-476DMR, all of which represent the clade-C sequence;

9. [5VN3](http://www.rcsb.org/pdb/explore/explore.do?structureID=5vn3) -- HIV-1 gp140 SOSIP trimer with bound sCD4 and FAB 17b; option to dock small-molecule CD4 mimetic BNM-III-170

10. [3G9R](http://www.rcsb/org/pdb/explore/explore.do?structureID=3g9r) -- HIV-1 gp41 membrane proximal external region (MPER) in an engineered coiled-coil trimer  

11. [5GGR](http://www.rcsb.org/pdb/explore/explore.do?structureID=5ggr) -- Nivolumab (``Opdivo'') Fab in complex with PD-1, with options for not including PD-1, and including either sucrose or D-mannitol exipients

12. [1L2Y](http://www.rcsb.org/pdb/explore/explore.do?structureID=1l2y) -- Trp cage miniprotein

13. SUCR -- a single solvated sucrose molecule extracted from 5o8l.pdb

14. [2JIU](http://www.rcsb.org/pdb/explore.do?structureID=2jiu) -- Human eGFR kinase, T790M mutant, ATPMg-bound

15. ALAD -- a single solvated alanine dipeptide with neutral ends (RESI ALAD in charmm36)

16. MTL -- a single solvated D-mannitol molecule extracted from 1m2w.pdb

17. B529 -- a single solvated BMS-529 molecule extracted from 5u7o.pdb, using CGenFF

18. [5U7O](http://www.rcsb.org/pdb/explore/explore.do?structureID=5u7o) -- HIV-1 gp140 SOSIP trimer with bound BMS-529 entry inhibitor, including glycans but not including Fabs.  There is also a version of this molecule with a representative BMS-derived DAVEI inhibitor model-built into the structure (5u7o-davei-l7).

19. [3PTB](http://www.rcsb.org/pdb/explore/explore.do?structureID=3ptb) -- Trypsin/benzamidine (latter parameterized using CGenFF)

20. [2EZN](http://www.rcsb.org/pdb/explore/explore.do?structureID=2ezn) -- Cyanovirin-N, with an option to create a CVN-(G4S)n-H6-MPER DAVEI

21. [2YHH](http://www.rcsb.org/pdb/explore/explore.do?structureID=2yhh) -- Microvirin, with options to create the MVN(Q81K/M83R)-(G4S)n-H6-[MPER|Trp3] DAVEI

22. BNM -- a single molecule of BNM-III-170 extracted from 5f4p.pdb, using CGenFF, in a water box

23. [5F4P](http://www.rcsb.org/pdb/explore/explore.do?structureID=5f4p) -- HIV-1 gp210 core with small-molecule CD4 mimic BNM-III-170 bound

24. [5JYN](http://www.rcsb.org/pdb/explore/explore.do?structureID=5jyn) -- HIV-1 gp41 transmembrane domain triple-helix embedded in a DMPC bilayer

25. [5VN8](http://www.rcsb.org/pdb/explore/explore.do?structureID=5vn8) -- HIV-1 gp140 SOSIP trimer in the b12-bound "open" conformation, with option to model-in MPER helices and the small-molecule entry-inhibitor BNM-III-170.

26. [2K7W](http://www.rcsb.org/pdb/explore/explore.do?structureID=2k7w) -- BAX proapoptotic protein with option to include bound BIM SAHB peptide.

27. [3CP1](http://www.rcsb.org/pdb/explore/explore.do?structureID=3cp1) -- HIV-1 gp41 NHR/CHR six-helix bundle, with option to grow in MPER and TM, and membrane-embedding.

28. EOH -- a single ethanol molecule extracted from 3TOD, using CGENFF

29. GXG -- a single GXG tripeptide, where X can be any residue, with 
an option to make a mixture of any molarity of GXG with an ethanol/water
cosolvent

31. [6VSB](http://www.rcsb.org/pdb/explore/explore.do?structureID=6vsb) -- Soluble, stabilized trimeric SARS-CoV-2 S spike protein complex, open

32. [6VXX](http://www.rcsb.org/pdb/explore/explore.do?structureID=6vxx) -- Soluble, stabilized trimeric SARS-CoV-2 S spike protein complex, closed

33. [6VYB](http://www.rcsb.org/pdb/explore/explore.do?structureID=6vyb) -- Soluble, stabilized trimeric SARS-CoV-2 S spike protein complex, open

34. [6M0J](http://www.rcsb.org/pdb/explore/explore/explore.do?structureID=6m0j) -- Complex of SARS-CoV-2 S spike receptor binding domain and ACE2

34. [6W41](http://www.rcsb.org/pdb/explore/explore/explore.do?structureID=6w41) -- Complex of SARS-CoV-2 S spike receptor binding domain and hAb CR3022

35. [6WAQ](http://www.rcsb.org/pdb/explore/explore/explore.do?structureID=6waq) -- Complex of SARS-CoV-1 S spike receptor binding domain and nanobody VHH-72

36. [4BYH](http://www.rcsb.org/pdb/explore/explore/explore.do?structureID=4byh) -- Sialylated IgG Fc

More to come...

## Acknowledgments

1. [VMD](http://www.ks.uiuc.edu/Research/vmd) and [NAMD](http://www.ks.uiuc.edu/Research/namd) were developed at the [Theoretical and Computational Biophysics Group at the NIH Center for Macromolecular Modeling and Bioinformatics at the University of Illinois at Urbana-Champaign](http://www.ks.uiuc.edu).  Please cite "W. Humphrey, A. Dalke, and K. Schulten.  VMD -- Visual Molecular Dynamics. Journal of Molecular Graphics, 1996;14:33-38" (VMD) and "J. C. Phillips, R. Braun, W. Wang, J. Gumbart, E. Tajkhorshid, E. Villa, C. Chipot, R. D. Skeel, L. Kale, and K. Schulten. Scalable molecular dynamics with NAMD. Journal of Computational Chemistry, 2005;26:1781-1802" (NAMD).

2. [Packmol](https://www.ime.unicamp.br/~martinez/packmol/userguide.shtml) is a product of Leandro Martinez in the Institute of Chemistry at the University of Campinas.  Please cite ``L. Martínez, R. Andrade, E. G. Birgin, J. M. Martínez. Packmol: A package for building initial configurations for molecular dynamics simulations. Journal of Computational Chemistry, 2009;30:2157-2164.'' 

3. The [CHARMM force field](http://mackerell.umaryland.edu/charmm_ff.shtml) is used in building all topology files in this repository.  Please cite: ``A. D. MacKerell, Jr., M. Feig, and C. L. Brooks, III. Extending the treatment of backbone energetics in protein force fields: limitations of gas-phase quantum mechanics in reproducing protein conformational distributions in molecular dynamics simulations. Journal of Computational Chemistry, 2004;25:1400-1415.''

4. [CGenFF](https://cgenff.paramchem.org/) is a web-based service for generating CHARMM General Forcefield parameters that are not already in the CHARMM36 force field.  Please cite: ``K. Vanommeslaeghe, E. Hatcher, C. Acharya, S. Kundu, S. Zhong, J. Shim, E. Darian, O. Guvench, P. Lopes, I. Vorobyov, A. D. MacKerell Jr., CHARMM General Force Field: A Force field for Drug-Like Molecules Compatible with the CHARMM All-Atom Additive Biological Force Field, Journal of Computational Chemistry, 2010;31,671-690.'' 

5. The [RCSB](https://rcsb.org) is funded by grant DBI-1338415 from the National Science Foundation, the National Institutes of Health, and the US Department of Energy.  Please cite:  ``H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne. The Protein Data Bank. Nucleic Acids Research, 2000;28:235-242. doi:10.1093/nar/28.1.235.''

6. All codes and data in this repository have been made possible with partial support from NIH through grants AI084117, AI093248, GM115249, GM056550, and GM100472, the National Science Foundation through grants DMR-1207389 and MCB-1330205, and the US Army through grants W911NF-12-2-0022, W911-NF-13-1-0046, W911NF-12-R-0011, and W911NF-17-2-0227.

2017-2020, Cameron F Abrams

cfa22@drexel.edu

