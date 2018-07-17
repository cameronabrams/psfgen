# 2yhh -- Microvirin ("CVN") with options to make MVN-(G4S)n-H6-[MPER|Trp3] DAVEI

## Files

This directory contains five files:
1. mkpsf_2yhh.tcl -- the VMD/psfgen script that creates the first vacuum psf/pdb pair.  
2. my_2yhh_vac.namd -- NAMD configuration file used to relax the "guessed" coordinates resulting from mkpsf_2yhh.tcl.
3. my_2yhh_solv.tcl -- VMD script that uses `solvate` and `autoionize` to generate a neutralized, solvated MD system using the coordinates from step 2 as input.  This also generates my_2yhh_colvars_op.inp, a colvars input file that defines collective variables that allow for center-of-mass restraint and an orientational restraint.  (Colvars restraints are NOT recommended for production MD.)
4. my_2yhh_solv_stageN.namd -- NAMD configuration file template to perform a minimization and short NpT MD of the raw solvated system.
5. do_test.sh -- a bash script that does everything

## Instructions

If you have cloned this repository, then make sure PSFGEN_BASEDIR resolves to the root directory of your local copy.  If you did not
clone the repository, you will have to figure this part out on your own.  It is also assumed below that CHARMRUN resolves to your local charmrun executable and NAMD2 resolves to your local NAMD2 executable.  For me, these are /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun and /home/cfa/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2.

To create the system, just invoke the do_test script in a clean directory:

```
$ mkdir 2yhh
$ cd 2yhh
$ $PSFGEN_BASEDIR/2yhh/do_test.sh
```

## DAVEI & Mutation (Q81K/M83R) Option

To generate a DAVEI molecule [a fusion of MVN, N repeats of (Gly4Ser), a His-tag, and either the 20-residue HIV-1 MPER segment (DKWASLWNWFEITEWLWYIK) or the9-residue "Trp3" segment (DKWASLWNW)], do this:

```
# mkdir my_davei
# cd my_davei
$ $PSFGEN_BASEDIR/2ezn/do_test.sh -psfgen_args -davei <N> [-trp3] [-mutate] [-seed #]
```

where `<N>` is replaced by the number of repeats of the (Gly4Ser) linker unit you want.  The optional `-trp3` flag, if present, directs the script to use the Trp3 sequence instead of the full MPER sequence, in building the DAVEI. The optional `-mutate` flag, if present, directs the script to generate the Q81K/M83R mutant protein.  The optional `-seed` flag allows the user to set the random-number generator seed.

2018, Cameron F Abrams
