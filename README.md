# psfgen -- a repository of advanced VMD/TcL and psfgen/TcL scripts for generating configurations and CHARMM36 topologies for use with NAMD

This repository contains psfgen scripts, TcL scripts for use in VMD, and some associated CHARMM36 topology and parameter files and NAMD config files used to generate _initial conditions_ for production MD simulations of various systems.  It should be helpful to anyone already familiar with using psfgen to build systems using CHARMM topologies who also needs access to some more advanced system-building capabilities than are available in most tutorials.  Some of the features this repository provides beyond what psfgen can easily do are the following:

* support for rudimentary loop model-building to include residues missing in the experimental PDB file but present in the crystallized protein sequence;
* support for glycans and non-covalently linked sugars;
* integration between solvation and initial MD simulation config file (easy transfer of box size);
* support for down-puckered prolines;
* (more to come)...

The _src_ directory contains files that can be "sourced" by psfgen scripts.  The _charmm_ directory contains some custom topologies and parameters derived from the latest charmm36 parameter set (July 2016).  Other directory names indicate the PDB entry for which the files contained therein are applicable.

## Requirements

1. NAMD v 2.12
  * Set environment variables CHARMRUN and NAMD2 to point to your system's charmrun and namd2; examples:
     * > CHARMRUN=${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun
     * > NAMD2=${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2
2. VMD v. 1.8.3
3. latest CHARMM36 topologies and parameters (toppar_c36_jul16.tgz) unpacked in ${HOME}/charmm/toppar
4. This repository cloned into ${HOME}/psfgen (or if somewhere else, point to it with the environment variable PSFGEN_BASEDIR)

## Currently supported PDB's:

1. [3TGQ](http://www.rcsb.org/pdb/explore/explore.do?structureId=3tgq) -- unliganded core monomeric HIV-1 gp120, with NAG's;

2. [4ZMJ](http://www.rcsb.org/pdb/explore/explore.do?structureId=4zmj) -- unliganded soluble trimeric HIV-1 Env gp140, with glycans;

3. [5FUU](http://www.rcsb.org/pdb/explore/explore.do?structureId=5fuu) -- soluble, cleaved JR-FL HIV-1 trimer, with glycans and with PGT151 antibodies removed;

4. [1HHP](http://www.rcsb.org/pdb/explore/explore.do?structureID=1hhp) -- apo, dimeric HIV-1 protease;

5. [2MB5](http://www.rcsb.org/pdb/explore/explore.do?structureID=2mb5) -- carbonmonoxymyoglobin;

6. [1F7A](http://www.rcsb.org/pdb/explore/explore.do?structureID=1f7a) -- HIV-1 protease with bound substrate peptide;

7. [1HIW](http://www.rcsb.org/pdb/explore/explore.do?structureID=1hiw) -- HIV-1 matrix (MA) trimer;

8. [4H8W](http://www.rcsb.org/pdb/explore/explore.do?structureID=4h8w) -- HIV-1 gp120 core clade A/E, with options to implement any number of the following mutations: (1) S375H (to bring back to clade A/E WT), (2) H61Y, (3) Q105H, (4) V108I, or (5) NIK474-476DMR, all of which represent the clade-C sequence;

9. [5VN3](http://www.rcsb.org/pdb/explore/explore.do?structureID=5vn3) -- HIV-1 gp140 SOSIP trimer with bound sCD4 and FAB 17b

10. [3G9R](http://www.rcsb/org/pdb/explore/explore.do?structureID=3g9r) -- HIV-1 gp41 membrane proximal external region (MPER) in an engineered coiled-coil trimer  

11. [5GGR](http://www.rcsb.org/pdb/explore/explore.do?structureID=5ggr) -- Nivolumab (``Opdivo'') Fab in complex with PD-1

12. [1L2Y](http://www.rcsb.org/pdb/explore/explore.do?structureID=1l2y) -- Trp cage miniprotein

13. [SUCR] -- a single solvated sucrose molecule extracted from 5o8l.pdb

More to come...

## Acknowledgments

1. VMD and NAMD are products of the [Theoretical and Computational Biophysics Group at the NIH Center for Macromolecular Modeling and Bioinformatics at the University of Illinois at Urbana-Champaign](http://www.ks.uiuc.edu)

2. All codes and data in this repository have been made possible with partial support from NIH through grants AI084117, AI093248, GM115249, GM056550, and GM100472, the National Science Foundation through grants DMR-1207389 and MCB-1330205, and the US Army through grants W911NF-12-2-0022, W911-NF-13-1-0046, and W911NF-12-R-0011.

2017-2018, Cameron F Abrams
cfa22@drexel.edu

