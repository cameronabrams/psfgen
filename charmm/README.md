# Some CHARMM-style topologies and parameter files

1. prod.top -- defines the "down-puckered" proline; CHARMM36 assumes all prolines are "up-puckered"
2. toppar_water_ions_namd_nonbfixes.str -- same as the toppar_water_ions_namd.str provided by A. MacKerrel except all NBFIXES are commented out.  You should copy this file to your ~/charmm/toppar/ directory, and only uses it in systems with proteins only.  I use this to prevent having to load every charmm36 parameter file, since it assumes the existence of atom types that are not in par_all36_prot.

Cameron F. Abrams, 2017
