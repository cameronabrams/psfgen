# VMD/psfgen script for generating psf/pdb pair for PDB 6vyb
# soluble, stabilized SARS-CoV-2 S trimer, symmetric, open
#
# cameron f abrams (c) 2020
# drexel university
# chemical and biological engineering
#
# check for base directory name variable;
# if not set, use default
if {![info exists PSFGEN_BASEDIR]} {
  # see if user set an environment variable
  if {[info exists env(PSFGEN_BASEDIR)]} {
      set PSFGEN_BASEDIR $env(PSFGEN_BASEDIR)
  } else {
      set PSFGEN_BASEDIR $env(HOME)/research/psfgen
  }
}

set seed 12345
set LOG_DCD 0
set logid -1
set P986K 0
set P987V 0
set CLEAVE 0
for { set a 0 } { $a < [llength $argv] } { incr a } {
  set arg [lindex $argv $a]
  if { $arg == "-seed" } {
    incr a
    set seed [lindex $argv $a]
  }
  if { $arg == "-log-dcd" } {
    set LOG_DCD 1
    incr a
    set log_dcd_file [lindex $argv $a]
  }
  if { $arg == "WT" } {
     set P986K 1
     set P987V 1
  }
  if { $arg == "CLEAVE" } {
     set CLEAVE 1
  }
}

expr srand($seed)

# load some custom TcL procedures to set coordinates correctly
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
set LOCALFILES {}

set DOMC 1

mol new 6vyb.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb_namd_cfa.rtf
topology $env(HOME)/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD

pdbalias residue NAG BGNA
pdbalias atom BGNA C7 C
pdbalias atom BGNA O7 O
pdbalias atom BGNA C8 CT

##### output of python3 parse_pdb_psfgen.py 6vyb.pdb below #####
set segs  { { A   27   69 }  { A   82  113 }  { A  116  143 }  { A  166  172 } 
            { A  186  242 }  { A  263  442 }  { A  448  470 }  { A  490  501 } 
            { A  503  620 }  { A  641  676 }  { A  690  811 }  { A  813  827 } 
            { A  855 1147 } 
            { B   27   66 }  { B   81  140 }  { B  164  172 }  { B  186  196 } 
            { B  200  211 }  { B  215  242 }  { B  263  454 }  { B  462  466 } 
            { B  491  515 }  { B  522  620 }  { B  641  676 }  { B  689  811 } 
            { B  813  827 }  { B  854 1147 } 
            { C   27   66 }  { C   81  143 }  { C  165  172 }  { C  186  242 } 
            { C  264  444 }  { C  448  454 }  { C  462  470 }  { C  491  620 } 
            { C  641  676 }  { C  690  811 }  { C  813  827 }  { C  856 1147 } 
            }
set loops { { A   70   81 }  { A  114  115 }  { A  144  165 }  { A  173  185 } 
            { A  243  262 }  { A  443  447 }  { A  471  489 }  { A  502  502 } 
            { A  621  640 }  { A  677  689 }  { A  812  812 }  { A  828  854 } 
           
            { B   67   80 }  { B  141  163 }  { B  173  185 }  { B  197  199 } 
            { B  212  214 }  { B  243  262 }  { B  455  461 }  { B  467  490 } 
            { B  516  521 }  { B  621  640 }  { B  677  688 }  { B  812  812 } 
            { B  828  853 } 
            { C   67   80 }  { C  144  164 }  { C  173  185 }  { C  243  263 } 
            { C  445  447 }  { C  455  461 }  { C  471  490 }  { C  621  640 } 
            { C  677  689 }  { C  812  812 }  { C  828  855 }  }
[atomselect top "chain A and protein and resid 27 to 69"] writepdb "A_27_to_69.pdb"
[atomselect top "chain A and protein and resid 82 to 113"] writepdb "A_82_to_113.pdb"
[atomselect top "chain A and protein and resid 116 to 143"] writepdb "A_116_to_143.pdb"
[atomselect top "chain A and protein and resid 166 to 172"] writepdb "A_166_to_172.pdb"
[atomselect top "chain A and protein and resid 186 to 242"] writepdb "A_186_to_242.pdb"
[atomselect top "chain A and protein and resid 263 to 442"] writepdb "A_263_to_442.pdb"
[atomselect top "chain A and protein and resid 448 to 470"] writepdb "A_448_to_470.pdb"
[atomselect top "chain A and protein and resid 490 to 501"] writepdb "A_490_to_501.pdb"
[atomselect top "chain A and protein and resid 503 to 620"] writepdb "A_503_to_620.pdb"
[atomselect top "chain A and protein and resid 641 to 676"] writepdb "A_641_to_676.pdb"
[atomselect top "chain A and protein and resid 690 to 811"] writepdb "A_690_to_811.pdb"
[atomselect top "chain A and protein and resid 813 to 827"] writepdb "A_813_to_827.pdb"
[atomselect top "chain A and protein and resid 855 to 1147"] writepdb "A_855_to_1147.pdb"
segment A {
   pdb A_27_to_69.pdb
   residue 70 VAL A
   residue 71 SER A
   residue 72 GLY A
   residue 73 THR A
   residue 74 ASN A
   residue 75 GLY A
   residue 76 THR A
   residue 77 LYS A
   residue 78 ARG A
   residue 79 PHE A
   residue 80 ASP A
   residue 81 ASN A
   pdb A_82_to_113.pdb
   residue 114 THR A
   residue 115 GLN A
   pdb A_116_to_143.pdb
   residue 144 TYR A
   residue 145 TYR A
   residue 146 HSE A
   residue 147 LYS A
   residue 148 ASN A
   residue 149 ASN A
   residue 150 LYS A
   residue 151 SER A
   residue 152 TRP A
   residue 153 MET A
   residue 154 GLU A
   residue 155 SER A
   residue 156 GLU A
   residue 157 PHE A
   residue 158 ARG A
   residue 159 VAL A
   residue 160 TYR A
   residue 161 SER A
   residue 162 SER A
   residue 163 ALA A
   residue 164 ASN A
   residue 165 ASN A
   pdb A_166_to_172.pdb
   residue 173 GLN A
   residue 174 PRO A
   residue 175 PHE A
   residue 176 LEU A
   residue 177 MET A
   residue 178 ASP A
   residue 179 LEU A
   residue 180 GLU A
   residue 181 GLY A
   residue 182 LYS A
   residue 183 GLN A
   residue 184 GLY A
   residue 185 ASN A
   pdb A_186_to_242.pdb
   residue 243 ALA A
   residue 244 LEU A
   residue 245 HSE A
   residue 246 ARG A
   residue 247 SER A
   residue 248 TYR A
   residue 249 LEU A
   residue 250 THR A
   residue 251 PRO A
   residue 252 GLY A
   residue 253 ASP A
   residue 254 SER A
   residue 255 SER A
   residue 256 SER A
   residue 257 GLY A
   residue 258 TRP A
   residue 259 THR A
   residue 260 ALA A
   residue 261 GLY A
   residue 262 ALA A
   pdb A_263_to_442.pdb
   residue 443 SER A
   residue 444 LYS A
   residue 445 VAL A
   residue 446 GLY A
   residue 447 GLY A
   pdb A_448_to_470.pdb
   residue 471 GLU A
   residue 472 ILE A
   residue 473 TYR A
   residue 474 GLN A
   residue 475 ALA A
   residue 476 GLY A
   residue 477 SER A
   residue 478 THR A
   residue 479 PRO A
   residue 480 CYS A
   residue 481 ASN A
   residue 482 GLY A
   residue 483 VAL A
   residue 484 GLU A
   residue 485 GLY A
   residue 486 PHE A
   residue 487 ASN A
   residue 488 CYS A
   residue 489 TYR A
   pdb A_490_to_501.pdb
   residue 502 GLY A
   pdb A_503_to_620.pdb
   residue 621 PRO A
   residue 622 VAL A
   residue 623 ALA A
   residue 624 ILE A
   residue 625 HSE A
   residue 626 ALA A
   residue 627 ASP A
   residue 628 GLN A
   residue 629 LEU A
   residue 630 THR A
   residue 631 PRO A
   residue 632 THR A
   residue 633 TRP A
   residue 634 ARG A
   residue 635 VAL A
   residue 636 TYR A
   residue 637 SER A
   residue 638 THR A
   residue 639 GLY A
   residue 640 SER A
   pdb A_641_to_676.pdb
   residue 677 GLN A
   residue 678 THR A
   residue 679 ASN A
   residue 680 SER A
   residue 681 PRO A
   residue 682 SER A
   residue 683 GLY A
   residue 684 ALA A
   residue 685 GLY A
   residue 686 SER A
   residue 687 VAL A
   residue 688 ALA A
   residue 689 SER A
   pdb A_690_to_811.pdb
   residue 812 PRO A
   pdb A_813_to_827.pdb
   residue 828 LEU A
   residue 829 ALA A
   residue 830 ASP A
   residue 831 ALA A
   residue 832 GLY A
   residue 833 PHE A
   residue 834 ILE A
   residue 835 LYS A
   residue 836 GLN A
   residue 837 TYR A
   residue 838 GLY A
   residue 839 ASP A
   residue 840 CYS A
   residue 841 LEU A
   residue 842 GLY A
   residue 843 ASP A
   residue 844 ILE A
   residue 845 ALA A
   residue 846 ALA A
   residue 847 ARG A
   residue 848 ASP A
   residue 849 LEU A
   residue 850 ILE A
   residue 851 CYS A
   residue 852 ALA A
   residue 853 GLN A
   residue 854 LYS A
   pdb A_855_to_1147.pdb
   if { $P986K == 1 } {
       mutate 986 LYS
   }
   if { $P987V == 1 } {
       mutate 987 VAL
   }
}
coordpdb A_27_to_69.pdb A
coordpdb A_82_to_113.pdb A
coordpdb A_116_to_143.pdb A
coordpdb A_166_to_172.pdb A
coordpdb A_186_to_242.pdb A
coordpdb A_263_to_442.pdb A
coordpdb A_448_to_470.pdb A
coordpdb A_490_to_501.pdb A
coordpdb A_503_to_620.pdb A
coordpdb A_641_to_676.pdb A
coordpdb A_690_to_811.pdb A
coordpdb A_813_to_827.pdb A
coordpdb A_855_to_1147.pdb A
coord A 70 N [cacoIn_nOut 69 A 0]
coord A 114 N [cacoIn_nOut 113 A 0]
coord A 144 N [cacoIn_nOut 143 A 0]
coord A 173 N [cacoIn_nOut 172 A 0]
coord A 243 N [cacoIn_nOut 242 A 0]
coord A 443 N [cacoIn_nOut 442 A 0]
coord A 471 N [cacoIn_nOut 470 A 0]
coord A 502 N [cacoIn_nOut 501 A 0]
coord A 621 N [cacoIn_nOut 620 A 0]
coord A 677 N [cacoIn_nOut 676 A 0]
coord A 812 N [cacoIn_nOut 811 A 0]
coord A 828 N [cacoIn_nOut 827 A 0]
[atomselect top "chain B and protein and resid 27 to 66"] writepdb "B_27_to_66.pdb"
[atomselect top "chain B and protein and resid 81 to 140"] writepdb "B_81_to_140.pdb"
[atomselect top "chain B and protein and resid 164 to 172"] writepdb "B_164_to_172.pdb"
[atomselect top "chain B and protein and resid 186 to 196"] writepdb "B_186_to_196.pdb"
[atomselect top "chain B and protein and resid 200 to 211"] writepdb "B_200_to_211.pdb"
[atomselect top "chain B and protein and resid 215 to 242"] writepdb "B_215_to_242.pdb"
[atomselect top "chain B and protein and resid 263 to 454"] writepdb "B_263_to_454.pdb"
[atomselect top "chain B and protein and resid 462 to 466"] writepdb "B_462_to_466.pdb"
[atomselect top "chain B and protein and resid 491 to 515"] writepdb "B_491_to_515.pdb"
[atomselect top "chain B and protein and resid 522 to 620"] writepdb "B_522_to_620.pdb"
[atomselect top "chain B and protein and resid 641 to 676"] writepdb "B_641_to_676.pdb"
[atomselect top "chain B and protein and resid 689 to 811"] writepdb "B_689_to_811.pdb"
[atomselect top "chain B and protein and resid 813 to 827"] writepdb "B_813_to_827.pdb"
[atomselect top "chain B and protein and resid 854 to 1147"] writepdb "B_854_to_1147.pdb"
segment B {
   pdb B_27_to_66.pdb
   residue 67 ALA B
   residue 68 ILE B
   residue 69 HSE B
   residue 70 VAL B
   residue 71 SER B
   residue 72 GLY B
   residue 73 THR B
   residue 74 ASN B
   residue 75 GLY B
   residue 76 THR B
   residue 77 LYS B
   residue 78 ARG B
   residue 79 PHE B
   residue 80 ASP B
   pdb B_81_to_140.pdb
   residue 141 LEU B
   residue 142 GLY B
   residue 143 VAL B
   residue 144 TYR B
   residue 145 TYR B
   residue 146 HSE B
   residue 147 LYS B
   residue 148 ASN B
   residue 149 ASN B
   residue 150 LYS B
   residue 151 SER B
   residue 152 TRP B
   residue 153 MET B
   residue 154 GLU B
   residue 155 SER B
   residue 156 GLU B
   residue 157 PHE B
   residue 158 ARG B
   residue 159 VAL B
   residue 160 TYR B
   residue 161 SER B
   residue 162 SER B
   residue 163 ALA B
   pdb B_164_to_172.pdb
   residue 173 GLN B
   residue 174 PRO B
   residue 175 PHE B
   residue 176 LEU B
   residue 177 MET B
   residue 178 ASP B
   residue 179 LEU B
   residue 180 GLU B
   residue 181 GLY B
   residue 182 LYS B
   residue 183 GLN B
   residue 184 GLY B
   residue 185 ASN B
   pdb B_186_to_196.pdb
   residue 197 ILE B
   residue 198 ASP B
   residue 199 GLY B
   pdb B_200_to_211.pdb
   residue 212 LEU B
   residue 213 VAL B
   residue 214 ARG B
   pdb B_215_to_242.pdb
   residue 243 ALA B
   residue 244 LEU B
   residue 245 HSE B
   residue 246 ARG B
   residue 247 SER B
   residue 248 TYR B
   residue 249 LEU B
   residue 250 THR B
   residue 251 PRO B
   residue 252 GLY B
   residue 253 ASP B
   residue 254 SER B
   residue 255 SER B
   residue 256 SER B
   residue 257 GLY B
   residue 258 TRP B
   residue 259 THR B
   residue 260 ALA B
   residue 261 GLY B
   residue 262 ALA B
   pdb B_263_to_454.pdb
   residue 455 LEU B
   residue 456 PHE B
   residue 457 ARG B
   residue 458 LYS B
   residue 459 SER B
   residue 460 ASN B
   residue 461 LEU B
   pdb B_462_to_466.pdb
   residue 467 ASP B
   residue 468 ILE B
   residue 469 SER B
   residue 470 THR B
   residue 471 GLU B
   residue 472 ILE B
   residue 473 TYR B
   residue 474 GLN B
   residue 475 ALA B
   residue 476 GLY B
   residue 477 SER B
   residue 478 THR B
   residue 479 PRO B
   residue 480 CYS B
   residue 481 ASN B
   residue 482 GLY B
   residue 483 VAL B
   residue 484 GLU B
   residue 485 GLY B
   residue 486 PHE B
   residue 487 ASN B
   residue 488 CYS B
   residue 489 TYR B
   residue 490 PHE B
   pdb B_491_to_515.pdb
   residue 516 GLU B
   residue 517 LEU B
   residue 518 LEU B
   residue 519 HSE B
   residue 520 ALA B
   residue 521 PRO B
   pdb B_522_to_620.pdb
   residue 621 PRO B
   residue 622 VAL B
   residue 623 ALA B
   residue 624 ILE B
   residue 625 HSE B
   residue 626 ALA B
   residue 627 ASP B
   residue 628 GLN B
   residue 629 LEU B
   residue 630 THR B
   residue 631 PRO B
   residue 632 THR B
   residue 633 TRP B
   residue 634 ARG B
   residue 635 VAL B
   residue 636 TYR B
   residue 637 SER B
   residue 638 THR B
   residue 639 GLY B
   residue 640 SER B
   pdb B_641_to_676.pdb
   residue 677 GLN B
   residue 678 THR B
   residue 679 ASN B
   residue 680 SER B
   residue 681 PRO B
   residue 682 SER B
   residue 683 GLY B
   residue 684 ALA B
   residue 685 GLY B
   residue 686 SER B
   residue 687 VAL B
   residue 688 ALA B
   pdb B_689_to_811.pdb
   residue 812 PRO B
   pdb B_813_to_827.pdb
   residue 828 LEU B
   residue 829 ALA B
   residue 830 ASP B
   residue 831 ALA B
   residue 832 GLY B
   residue 833 PHE B
   residue 834 ILE B
   residue 835 LYS B
   residue 836 GLN B
   residue 837 TYR B
   residue 838 GLY B
   residue 839 ASP B
   residue 840 CYS B
   residue 841 LEU B
   residue 842 GLY B
   residue 843 ASP B
   residue 844 ILE B
   residue 845 ALA B
   residue 846 ALA B
   residue 847 ARG B
   residue 848 ASP B
   residue 849 LEU B
   residue 850 ILE B
   residue 851 CYS B
   residue 852 ALA B
   residue 853 GLN B
   pdb B_854_to_1147.pdb
   if { $P986K == 1 } {
       mutate 986 LYS
   }
   if { $P987V == 1 } {
       mutate 987 VAL
   }
}
coordpdb B_27_to_66.pdb B
coordpdb B_81_to_140.pdb B
coordpdb B_164_to_172.pdb B
coordpdb B_186_to_196.pdb B
coordpdb B_200_to_211.pdb B
coordpdb B_215_to_242.pdb B
coordpdb B_263_to_454.pdb B
coordpdb B_462_to_466.pdb B
coordpdb B_491_to_515.pdb B
coordpdb B_522_to_620.pdb B
coordpdb B_641_to_676.pdb B
coordpdb B_689_to_811.pdb B
coordpdb B_813_to_827.pdb B
coordpdb B_854_to_1147.pdb B
coord B 67 N [cacoIn_nOut 66 B 0]
coord B 141 N [cacoIn_nOut 140 B 0]
coord B 173 N [cacoIn_nOut 172 B 0]
coord B 197 N [cacoIn_nOut 196 B 0]
coord B 212 N [cacoIn_nOut 211 B 0]
coord B 243 N [cacoIn_nOut 242 B 0]
coord B 455 N [cacoIn_nOut 454 B 0]
coord B 467 N [cacoIn_nOut 466 B 0]
coord B 516 N [cacoIn_nOut 515 B 0]
coord B 621 N [cacoIn_nOut 620 B 0]
coord B 677 N [cacoIn_nOut 676 B 0]
coord B 812 N [cacoIn_nOut 811 B 0]
coord B 828 N [cacoIn_nOut 827 B 0]
[atomselect top "chain C and protein and resid 27 to 66"] writepdb "C_27_to_66.pdb"
[atomselect top "chain C and protein and resid 81 to 143"] writepdb "C_81_to_143.pdb"
[atomselect top "chain C and protein and resid 165 to 172"] writepdb "C_165_to_172.pdb"
[atomselect top "chain C and protein and resid 186 to 242"] writepdb "C_186_to_242.pdb"
[atomselect top "chain C and protein and resid 264 to 444"] writepdb "C_264_to_444.pdb"
[atomselect top "chain C and protein and resid 448 to 454"] writepdb "C_448_to_454.pdb"
[atomselect top "chain C and protein and resid 462 to 470"] writepdb "C_462_to_470.pdb"
[atomselect top "chain C and protein and resid 491 to 620"] writepdb "C_491_to_620.pdb"
[atomselect top "chain C and protein and resid 641 to 676"] writepdb "C_641_to_676.pdb"
[atomselect top "chain C and protein and resid 690 to 811"] writepdb "C_690_to_811.pdb"
[atomselect top "chain C and protein and resid 813 to 827"] writepdb "C_813_to_827.pdb"
[atomselect top "chain C and protein and resid 856 to 1147"] writepdb "C_856_to_1147.pdb"
segment C {
   pdb C_27_to_66.pdb
   residue 67 ALA C
   residue 68 ILE C
   residue 69 HSE C
   residue 70 VAL C
   residue 71 SER C
   residue 72 GLY C
   residue 73 THR C
   residue 74 ASN C
   residue 75 GLY C
   residue 76 THR C
   residue 77 LYS C
   residue 78 ARG C
   residue 79 PHE C
   residue 80 ASP C
   pdb C_81_to_143.pdb
   residue 144 TYR C
   residue 145 TYR C
   residue 146 HSE C
   residue 147 LYS C
   residue 148 ASN C
   residue 149 ASN C
   residue 150 LYS C
   residue 151 SER C
   residue 152 TRP C
   residue 153 MET C
   residue 154 GLU C
   residue 155 SER C
   residue 156 GLU C
   residue 157 PHE C
   residue 158 ARG C
   residue 159 VAL C
   residue 160 TYR C
   residue 161 SER C
   residue 162 SER C
   residue 163 ALA C
   residue 164 ASN C
   pdb C_165_to_172.pdb
   residue 173 GLN C
   residue 174 PRO C
   residue 175 PHE C
   residue 176 LEU C
   residue 177 MET C
   residue 178 ASP C
   residue 179 LEU C
   residue 180 GLU C
   residue 181 GLY C
   residue 182 LYS C
   residue 183 GLN C
   residue 184 GLY C
   residue 185 ASN C
   pdb C_186_to_242.pdb
   residue 243 ALA C
   residue 244 LEU C
   residue 245 HSE C
   residue 246 ARG C
   residue 247 SER C
   residue 248 TYR C
   residue 249 LEU C
   residue 250 THR C
   residue 251 PRO C
   residue 252 GLY C
   residue 253 ASP C
   residue 254 SER C
   residue 255 SER C
   residue 256 SER C
   residue 257 GLY C
   residue 258 TRP C
   residue 259 THR C
   residue 260 ALA C
   residue 261 GLY C
   residue 262 ALA C
   residue 263 ALA C
   pdb C_264_to_444.pdb
   residue 445 VAL C
   residue 446 GLY C
   residue 447 GLY C
   pdb C_448_to_454.pdb
   residue 455 LEU C
   residue 456 PHE C
   residue 457 ARG C
   residue 458 LYS C
   residue 459 SER C
   residue 460 ASN C
   residue 461 LEU C
   pdb C_462_to_470.pdb
   residue 471 GLU C
   residue 472 ILE C
   residue 473 TYR C
   residue 474 GLN C
   residue 475 ALA C
   residue 476 GLY C
   residue 477 SER C
   residue 478 THR C
   residue 479 PRO C
   residue 480 CYS C
   residue 481 ASN C
   residue 482 GLY C
   residue 483 VAL C
   residue 484 GLU C
   residue 485 GLY C
   residue 486 PHE C
   residue 487 ASN C
   residue 488 CYS C
   residue 489 TYR C
   residue 490 PHE C
   pdb C_491_to_620.pdb
   residue 621 PRO C
   residue 622 VAL C
   residue 623 ALA C
   residue 624 ILE C
   residue 625 HSE C
   residue 626 ALA C
   residue 627 ASP C
   residue 628 GLN C
   residue 629 LEU C
   residue 630 THR C
   residue 631 PRO C
   residue 632 THR C
   residue 633 TRP C
   residue 634 ARG C
   residue 635 VAL C
   residue 636 TYR C
   residue 637 SER C
   residue 638 THR C
   residue 639 GLY C
   residue 640 SER C
   pdb C_641_to_676.pdb
   residue 677 GLN C
   residue 678 THR C
   residue 679 ASN C
   residue 680 SER C
   residue 681 PRO C
   residue 682 SER C
   residue 683 GLY C
   residue 684 ALA C
   residue 685 GLY C
   residue 686 SER C
   residue 687 VAL C
   residue 688 ALA C
   residue 689 SER C
   pdb C_690_to_811.pdb
   residue 812 PRO C
   pdb C_813_to_827.pdb
   residue 828 LEU C
   residue 829 ALA C
   residue 830 ASP C
   residue 831 ALA C
   residue 832 GLY C
   residue 833 PHE C
   residue 834 ILE C
   residue 835 LYS C
   residue 836 GLN C
   residue 837 TYR C
   residue 838 GLY C
   residue 839 ASP C
   residue 840 CYS C
   residue 841 LEU C
   residue 842 GLY C
   residue 843 ASP C
   residue 844 ILE C
   residue 845 ALA C
   residue 846 ALA C
   residue 847 ARG C
   residue 848 ASP C
   residue 849 LEU C
   residue 850 ILE C
   residue 851 CYS C
   residue 852 ALA C
   residue 853 GLN C
   residue 854 LYS C
   residue 855 PHE C
   pdb C_856_to_1147.pdb
   if { $P986K == 1 } {
       mutate 986 LYS
   }
   if { $P987V == 1 } {
       mutate 987 VAL
   }
}
coordpdb C_27_to_66.pdb C
coordpdb C_81_to_143.pdb C
coordpdb C_165_to_172.pdb C
coordpdb C_186_to_242.pdb C
coordpdb C_264_to_444.pdb C
coordpdb C_448_to_454.pdb C
coordpdb C_462_to_470.pdb C
coordpdb C_491_to_620.pdb C
coordpdb C_641_to_676.pdb C
coordpdb C_690_to_811.pdb C
coordpdb C_813_to_827.pdb C
coordpdb C_856_to_1147.pdb C
coord C 67 N [cacoIn_nOut 66 C 0]
coord C 144 N [cacoIn_nOut 143 C 0]
coord C 173 N [cacoIn_nOut 172 C 0]
coord C 243 N [cacoIn_nOut 242 C 0]
coord C 445 N [cacoIn_nOut 444 C 0]
coord C 455 N [cacoIn_nOut 454 C 0]
coord C 471 N [cacoIn_nOut 470 C 0]
coord C 621 N [cacoIn_nOut 620 C 0]
coord C 677 N [cacoIn_nOut 676 C 0]
coord C 812 N [cacoIn_nOut 811 C 0]
coord C 828 N [cacoIn_nOut 827 C 0]
set myseg [atomselect top "chain A and resid 1301 to 1320"]
$myseg set resname BGNA
$myseg writepdb AS_1301_to_1320.pdb
segment AS {
   pdb AS_1301_to_1320.pdb
}
coordpdb AS_1301_to_1320.pdb AS
set myseg [atomselect top "chain B and resid 1301 to 1319"]
$myseg set resname BGNA
$myseg writepdb BS_1301_to_1319.pdb
segment BS {
   pdb BS_1301_to_1319.pdb
}
coordpdb BS_1301_to_1319.pdb BS
set myseg [atomselect top "chain C and resid 1301 to 1320"]
$myseg set resname BGNA
$myseg writepdb CS_1301_to_1320.pdb
segment CS {
   pdb CS_1301_to_1320.pdb
}
coordpdb CS_1301_to_1320.pdb CS
patch DISU A:131 A:166
patch DISU A:291 A:301
patch DISU A:336 A:361
patch DISU A:379 A:432
patch DISU A:391 A:525
patch DISU A:538 A:590
patch DISU A:617 A:649
patch DISU A:662 A:671
patch DISU A:738 A:760
patch DISU A:743 A:749
patch DISU A:1032 A:1043
patch DISU A:1082 A:1126
patch DISU B:131 B:166
patch DISU B:291 B:301
patch DISU B:379 B:432
patch DISU B:391 B:525
patch DISU B:538 B:590
patch DISU B:617 B:649
patch DISU B:662 B:671
patch DISU B:738 B:760
patch DISU B:743 B:749
patch DISU B:1032 B:1043
patch DISU B:1082 B:1126
patch DISU C:131 C:166
patch DISU C:291 C:301
patch DISU C:336 C:361
patch DISU C:379 C:432
patch DISU C:391 C:525
patch DISU C:538 C:590
patch DISU C:617 C:649
patch DISU C:662 C:671
patch DISU C:738 C:760
patch DISU C:743 C:749
patch DISU C:1032 C:1043
patch DISU C:1082 C:1126
patch NGLB A:61 AS:1301
patch NGLB A:122 AS:1302
patch NGLB A:234 AS:1303
patch NGLB A:282 AS:1305
patch NGLB A:331 AS:1306
patch NGLB A:343 AS:1307
patch NGLB A:603 AS:1308
patch NGLB A:616 AS:1309
patch NGLB A:657 AS:1310
patch NGLB A:709 AS:1311
patch NGLB A:717 AS:1312
patch NGLB A:801 AS:1314
patch NGLB A:1074 AS:1316
patch NGLB A:1098 AS:1317
patch NGLB A:1134 AS:1319
patch NGLB B:61 BS:1301
patch NGLB B:122 BS:1302
patch NGLB B:165 BS:1319
patch NGLB B:234 BS:1303
patch NGLB B:282 BS:1304
patch NGLB B:331 BS:1305
patch NGLB B:343 BS:1306
patch NGLB B:603 BS:1307
patch NGLB B:616 BS:1308
patch NGLB B:657 BS:1309
patch NGLB B:709 BS:1310
patch NGLB B:717 BS:1311
patch NGLB B:801 BS:1312
patch NGLB B:1074 BS:1314
patch NGLB B:1098 BS:1315
patch NGLB B:1134 BS:1317
patch NGLB C:61 CS:1301
patch NGLB C:122 CS:1302
patch NGLB C:165 CS:1320
patch NGLB C:234 CS:1303
patch NGLB C:282 CS:1304
patch NGLB C:331 CS:1305
patch NGLB C:343 CS:1306
patch NGLB C:603 CS:1307
patch NGLB C:616 CS:1308
patch NGLB C:657 CS:1309
patch NGLB C:709 CS:1310
patch NGLB C:717 CS:1311
patch NGLB C:801 CS:1313
patch NGLB C:1074 CS:1315
patch NGLB C:1098 CS:1316
patch NGLB C:1134 CS:1318
patch 14bb AS:1303 AS:1304
patch 14bb AS:1312 AS:1313
patch 14bb AS:1314 AS:1315
patch 14bb AS:1317 AS:1318
patch 14bb AS:1319 AS:1320
patch 14bb BS:1312 BS:1313
patch 14bb BS:1315 BS:1316
patch 14bb BS:1317 BS:1318
patch 14bb CS:1311 CS:1312
patch 14bb CS:1313 CS:1314
patch 14bb CS:1316 CS:1317
patch 14bb CS:1318 CS:1319
##### output of python3 parse_pdb_psfgen.py 6vyb.pdb above #####

guesscoord

regenerate angles dihedrals

writepsf "my_6vyb.psf"
writepdb "unrelaxed.pdb"

#lappend LOCALFILES unrelaxed.pdb

mol delete top
mol new my_6vyb.psf
mol addfile unrelaxed.pdb
set molid [molinfo top get id]
set or [measure center [atomselect top "all"] weight mass]
set a [atomselect top all]
$a moveby [vecscale -1 $or]
if { $LOG_DCD != "0" } {
   mol new my_6vyb.psf
   mol addfile unrelaxed.pdb
   set logid [molinfo top get id]
   mol top $molid
   log_addframe ${molid} ${logid}
}
if { 0 } {
set ca [measure center [atomselect top "protein and chain A C E"] weight mass]
set cb [measure center [atomselect top "protein and chain B D F"] weight mass]   
set pi 3.1415928
set dv [vecsub $ca $cb]
set d [veclength $dv]
set cp [expr [lindex $dv 0]/$d]
set sp [expr [lindex $dv 1]/$d]
set p [expr acos($cp)]
if {[expr $sp < 0.0]} {
  set p [expr 2*$pi-$p]
}
set ct [expr [lindex $dv 2]/$d]
set t [expr acos($ct)]
$a move [transaxis z [expr -1 * $p] rad]
$a move [transaxis y [expr -1 * $t] rad]
$a writepdb "unrelaxed2.pdb"
lappend LOCALFILES "unrelaxed2.pdb"
}

set nc 1000
set rcut 3.0
set temperature 3.0
set k 10.0
set r0 1.5
set bg [atomselect ${molid} "noh"]
set loopindex 0
set nloops [llength $loops]
foreach l $loops {
  set chain [lindex $l 0]
  puts "Relaxing loop $loopindex ($l) out of $nloops"
  set residueList [[atomselect ${molid} "chain $chain and resid [lindex $l 1] to [lindex $l 2] and name CA"] get residue]
  do_loop_mc ${residueList} ${chain} ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999] $logid
  set loopindex [expr $loopindex + 1]
}
$a writepdb "my_6vyb_mcOut.pdb"

if { $LOG_DCD != "0" } {
   set loga [atomselect $logid all]
   animate write dcd $log_dcd_file waitfor all sel $loga $logid
}

# make a pdb file that fixes all heavy atoms in the original
# crystal structure -- all added atoms are set as unfixed
# for a minimization
mol delete top
mol new my_6vyb.psf
mol addfile unrelaxed.pdb
set a [atomselect top all]
$a set beta 0
foreach s $segs {
  [atomselect top "chain [lindex $s 0] and resid [lindex $s 1] to [lindex $s 2]"] set beta 1
}
[atomselect top "not noh"] set beta 0

$a writepdb "my_6vyb_fix.pdb"

if { $CLEAVE == 1 } {
   mol delete top
   psfcontext reset
   puts "MKPSF> Sourcing cleave.tcl..."
   source $PSFGEN_BASEDIR/6vyb/cleave.tcl
}

# clean up
foreach f $LOCALFILES { 
  exec rm $f
}

quit

