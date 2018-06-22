# Some CHARMM-style topologies and parameter files

1. prod.top -- defines the "down-puckered" proline; CHARMM36 assumes all prolines are "up-puckered"
2. toppar_water_ions_namd_nonbfixes.str -- same as the toppar_water_ions_namd.str provided by A. MacKerrel except all NBFIXES are commented out.  You should copy this file to your ~/charmm/toppar/ directory, and only use it in systems with proteins only.  I use this to prevent having to load every charmm36 parameter file, since it assumes the existence of atom types that are not in par_all36_prot.
3. top_all36_carb_namd_cfa.rtf -- same as top_all36_carb.rtf except the RESIdue BGLCNA is changed to BGNA since evidently NAMD can't handle residues longer than four characters.  You should just copy this to your ~/charmm/toppar/ directory.
4. bms529.str -- Stream file from CGenFF for BMS-529 entry inhibitor
5. aeg -- Stream file for the AEG compound, a derivative of bms529 with a furanyl in place of the azole
6. dls1 and dls2 -- stream files for linker segments in the DAVEI compounds

Cameron F. Abrams, cfa22@drexel.edu, 2018
