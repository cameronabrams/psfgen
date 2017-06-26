# psfgen
## TcL, VMD and psfgen scripts for generating MD systems for use with NAMD

This repository contains some psfgen scripts, TcL scripts for use in VMD, and some associated CHARMM topology and parameter files and NAMD config files used to generate _initial conditions_ for production MD simulations of various protein systems.  It should be helpful to anyone already familiar with using psfgen to build systems.  Some of the features this repository provides beyond what psfgen can easily do:

* support for rudimentary loop model-building to include residues missing in the experimental PDB file but present in the crystallized protein sequence;  
* support for including glycans;
* integration between solvation and initial MD simulation config file (easy transfer of box size);
* support for down-puckered prolines;
* (more to come)...

The tcl directory contains files than can be "sourced" by psfgen scripts.  The charmm directory contains some custom topologies and parameters derived from the latest charmm36 parameter set (July 2016).  Other directory names indicate the PDB entry for which the files contained therein are applicable.

## Requirements

1. NAMD v 2.12
2. VMD v. 1.8.3
3. latest CHARMM36 topologies and parameters in $HOME/charmm/toppar

## Currently supported PDB's:

1. 3TGQ -- unliganded core monomeric HIV-1 gp120 [http://www.rcsb.org/pdb/explore/explore.do?structureId=3tgq]

2017, Cameron F Abrams
