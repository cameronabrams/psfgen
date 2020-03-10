# VMD/psfgen script for generating psf/pdb pair for PDB 6vsb
# soluble, stabilized SARS-CoV-2 S trimer
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
}

expr srand($seed)

# load some custom TcL procedures to set coordinates correctly
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
set LOCALFILES {}

set DOMC 1

mol new 6vsb.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb_namd_cfa.rtf
topology $env(HOME)/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD

#pdbalias residue NAG BGLC
#pdbalias atom BGLC C7 C
#pdbalias atom BGLC O7 O
#pdbalias atom BGLC C8 CT

##### output of python3 parse_pdb_psfgen.py 6vsb.pdb below ###
set segs  { { A   27   66 }  { A   79   95 }  { A   99  142 }  { A  156  176 } 
            { A  187  246 }  { A  261  328 }  { A  335  443 }  { A  449  454 } 
            { A  491  500 }  { A  503  620 }  { A  640  672 }  { A  687  811 } 
            { A  815  828 }  { A  853 1146 } 
            { B   27   66 }  { B   81  141 }  { B  155  176 }  { B  187  209 } 
            { B  217  242 }  { B  263  443 }  { B  449  454 }  { B  460  473 } 
            { B  487  500 }  { B  503  620 }  { B  638  672 }  { B  687  811 } 
            { B  815  828 }  { B  853 1146 } 
            { C   27   66 }  { C   80   95 }  { C   99  140 }  { C  157  176 } 
            { C  187  245 }  { C  261  443 }  { C  449  454 }  { C  460  471 } 
            { C  487  498 }  { C  503  620 }  { C  641  672 }  { C  687  811 } 
            { C  815  828 }  { C  853 1146 }  }
set loops { { A   67   78 }  { A   96   98 }  { A  143  155 }  { A  177  186 } 
            { A  247  260 }  { A  329  334 }  { A  444  448 }  { A  455  490 } 
            { A  501  502 }  { A  621  639 }  { A  673  686 }  { A  812  814 } 
            { A  829  852 } 
            { B   67   80 }  { B  142  154 }  { B  177  186 }  { B  210  216 } 
            { B  243  262 }  { B  444  448 }  { B  455  459 }  { B  474  486 } 
            { B  501  502 }  { B  621  637 }  { B  673  686 }  { B  812  814 } 
            { B  829  852 } 
            { C   67   79 }  { C   96   98 }  { C  141  156 }  { C  177  186 } 
            { C  246  260 }  { C  444  448 }  { C  455  459 }  { C  472  486 } 
            { C  499  502 }  { C  621  640 }  { C  673  686 }  { C  812  814 } 
            { C  829  852 }  }
[atomselect top "chain A and protein and resid 27 to 66"] writepdb "A_27_to_66.pdb"
[atomselect top "chain A and protein and resid 79 to 95"] writepdb "A_79_to_95.pdb"
[atomselect top "chain A and protein and resid 99 to 142"] writepdb "A_99_to_142.pdb"
[atomselect top "chain A and protein and resid 156 to 176"] writepdb "A_156_to_176.pdb"
[atomselect top "chain A and protein and resid 187 to 246"] writepdb "A_187_to_246.pdb"
[atomselect top "chain A and protein and resid 261 to 328"] writepdb "A_261_to_328.pdb"
[atomselect top "chain A and protein and resid 335 to 443"] writepdb "A_335_to_443.pdb"
[atomselect top "chain A and protein and resid 449 to 454"] writepdb "A_449_to_454.pdb"
[atomselect top "chain A and protein and resid 491 to 500"] writepdb "A_491_to_500.pdb"
[atomselect top "chain A and protein and resid 503 to 620"] writepdb "A_503_to_620.pdb"
[atomselect top "chain A and protein and resid 640 to 672"] writepdb "A_640_to_672.pdb"
[atomselect top "chain A and protein and resid 687 to 811"] writepdb "A_687_to_811.pdb"
[atomselect top "chain A and protein and resid 815 to 828"] writepdb "A_815_to_828.pdb"
[atomselect top "chain A and protein and resid 853 to 1146"] writepdb "A_853_to_1146.pdb"
segment A {
   pdb A_27_to_66.pdb
   residue 67 ALA A
   residue 68 ILE A
   residue 69 HSE A
   residue 70 VAL A
   residue 71 SER A
   residue 72 GLY A
   residue 73 THR A
   residue 74 ASN A
   residue 75 GLY A
   residue 76 THR A
   residue 77 LYS A
   residue 78 ARG A
   pdb A_79_to_95.pdb
   residue 96 GLU A
   residue 97 LYS A
   residue 98 SER A
   pdb A_99_to_142.pdb
   residue 143 VAL A
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
   pdb A_156_to_176.pdb
   residue 177 MET A
   residue 178 ASP A
   residue 179 LEU A
   residue 180 GLU A
   residue 181 GLY A
   residue 182 LYS A
   residue 183 GLN A
   residue 184 GLY A
   residue 185 ASN A
   residue 186 PHE A
   pdb A_187_to_246.pdb
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
   pdb A_261_to_328.pdb
   residue 329 PHE A
   residue 330 PRO A
   residue 331 ASN A
   residue 332 ILE A
   residue 333 THR A
   residue 334 ASN A
   pdb A_335_to_443.pdb
   residue 444 LYS A
   residue 445 VAL A
   residue 446 GLY A
   residue 447 GLY A
   residue 448 ASN A
   pdb A_449_to_454.pdb
   residue 455 LEU A
   residue 456 PHE A
   residue 457 ARG A
   residue 458 LYS A
   residue 459 SER A
   residue 460 ASN A
   residue 461 LEU A
   residue 462 LYS A
   residue 463 PRO A
   residue 464 PHE A
   residue 465 GLU A
   residue 466 ARG A
   residue 467 ASP A
   residue 468 ILE A
   residue 469 SER A
   residue 470 THR A
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
   residue 490 PHE A
   pdb A_491_to_500.pdb
   residue 501 ASN A
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
   pdb A_640_to_672.pdb
   residue 673 SER A
   residue 674 TYR A
   residue 675 GLN A
   residue 676 THR A
   residue 677 GLN A
   residue 678 THR A
   residue 679 ASN A
   residue 680 SER A
   residue 681 PRO A
   residue 682 GLY A
   residue 683 SER A
   residue 684 ALA A
   residue 685 SER A
   residue 686 SER A
   pdb A_687_to_811.pdb
   residue 812 PRO A
   residue 813 SER A
   residue 814 LYS A
   pdb A_815_to_828.pdb
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
   pdb A_853_to_1146.pdb
}
coordpdb A_27_to_66.pdb A
coordpdb A_79_to_95.pdb A
coordpdb A_99_to_142.pdb A
coordpdb A_156_to_176.pdb A
coordpdb A_187_to_246.pdb A
coordpdb A_261_to_328.pdb A
coordpdb A_335_to_443.pdb A
coordpdb A_449_to_454.pdb A
coordpdb A_491_to_500.pdb A
coordpdb A_503_to_620.pdb A
coordpdb A_640_to_672.pdb A
coordpdb A_687_to_811.pdb A
coordpdb A_815_to_828.pdb A
coordpdb A_853_to_1146.pdb A
coord A 67 N [cacoIn_nOut 66 A 0]
coord A 96 N [cacoIn_nOut 95 A 0]
coord A 143 N [cacoIn_nOut 142 A 0]
coord A 177 N [cacoIn_nOut 176 A 0]
coord A 247 N [cacoIn_nOut 246 A 0]
coord A 329 N [cacoIn_nOut 328 A 0]
coord A 444 N [cacoIn_nOut 443 A 0]
coord A 455 N [cacoIn_nOut 454 A 0]
coord A 501 N [cacoIn_nOut 500 A 0]
coord A 621 N [cacoIn_nOut 620 A 0]
coord A 673 N [cacoIn_nOut 672 A 0]
coord A 812 N [cacoIn_nOut 811 A 0]
coord A 829 N [cacoIn_nOut 828 A 0]
[atomselect top "chain B and protein and resid 27 to 66"] writepdb "B_27_to_66.pdb"
[atomselect top "chain B and protein and resid 81 to 141"] writepdb "B_81_to_141.pdb"
[atomselect top "chain B and protein and resid 155 to 176"] writepdb "B_155_to_176.pdb"
[atomselect top "chain B and protein and resid 187 to 209"] writepdb "B_187_to_209.pdb"
[atomselect top "chain B and protein and resid 217 to 242"] writepdb "B_217_to_242.pdb"
[atomselect top "chain B and protein and resid 263 to 443"] writepdb "B_263_to_443.pdb"
[atomselect top "chain B and protein and resid 449 to 454"] writepdb "B_449_to_454.pdb"
[atomselect top "chain B and protein and resid 460 to 473"] writepdb "B_460_to_473.pdb"
[atomselect top "chain B and protein and resid 487 to 500"] writepdb "B_487_to_500.pdb"
[atomselect top "chain B and protein and resid 503 to 620"] writepdb "B_503_to_620.pdb"
[atomselect top "chain B and protein and resid 638 to 672"] writepdb "B_638_to_672.pdb"
[atomselect top "chain B and protein and resid 687 to 811"] writepdb "B_687_to_811.pdb"
[atomselect top "chain B and protein and resid 815 to 828"] writepdb "B_815_to_828.pdb"
[atomselect top "chain B and protein and resid 853 to 1146"] writepdb "B_853_to_1146.pdb"
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
   pdb B_81_to_141.pdb
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
   pdb B_155_to_176.pdb
   residue 177 MET B
   residue 178 ASP B
   residue 179 LEU B
   residue 180 GLU B
   residue 181 GLY B
   residue 182 LYS B
   residue 183 GLN B
   residue 184 GLY B
   residue 185 ASN B
   residue 186 PHE B
   pdb B_187_to_209.pdb
   residue 210 ILE B
   residue 211 ASN B
   residue 212 LEU B
   residue 213 VAL B
   residue 214 ARG B
   residue 215 ASP B
   residue 216 LEU B
   pdb B_217_to_242.pdb
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
   pdb B_263_to_443.pdb
   residue 444 LYS B
   residue 445 VAL B
   residue 446 GLY B
   residue 447 GLY B
   residue 448 ASN B
   pdb B_449_to_454.pdb
   residue 455 LEU B
   residue 456 PHE B
   residue 457 ARG B
   residue 458 LYS B
   residue 459 SER B
   pdb B_460_to_473.pdb
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
   pdb B_487_to_500.pdb
   residue 501 ASN B
   residue 502 GLY B
   pdb B_503_to_620.pdb
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
   pdb B_638_to_672.pdb
   residue 673 SER B
   residue 674 TYR B
   residue 675 GLN B
   residue 676 THR B
   residue 677 GLN B
   residue 678 THR B
   residue 679 ASN B
   residue 680 SER B
   residue 681 PRO B
   residue 682 GLY B
   residue 683 SER B
   residue 684 ALA B
   residue 685 SER B
   residue 686 SER B
   pdb B_687_to_811.pdb
   residue 812 PRO B
   residue 813 SER B
   residue 814 LYS B
   pdb B_815_to_828.pdb
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
   pdb B_853_to_1146.pdb
}
coordpdb B_27_to_66.pdb B
coordpdb B_81_to_141.pdb B
coordpdb B_155_to_176.pdb B
coordpdb B_187_to_209.pdb B
coordpdb B_217_to_242.pdb B
coordpdb B_263_to_443.pdb B
coordpdb B_449_to_454.pdb B
coordpdb B_460_to_473.pdb B
coordpdb B_487_to_500.pdb B
coordpdb B_503_to_620.pdb B
coordpdb B_638_to_672.pdb B
coordpdb B_687_to_811.pdb B
coordpdb B_815_to_828.pdb B
coordpdb B_853_to_1146.pdb B
coord B 67 N [cacoIn_nOut 66 B 0]
coord B 142 N [cacoIn_nOut 141 B 0]
coord B 177 N [cacoIn_nOut 176 B 0]
coord B 210 N [cacoIn_nOut 209 B 0]
coord B 243 N [cacoIn_nOut 242 B 0]
coord B 444 N [cacoIn_nOut 443 B 0]
coord B 455 N [cacoIn_nOut 454 B 0]
coord B 474 N [cacoIn_nOut 473 B 0]
coord B 501 N [cacoIn_nOut 500 B 0]
coord B 621 N [cacoIn_nOut 620 B 0]
coord B 673 N [cacoIn_nOut 672 B 0]
coord B 812 N [cacoIn_nOut 811 B 0]
coord B 829 N [cacoIn_nOut 828 B 0]
[atomselect top "chain C and protein and resid 27 to 66"] writepdb "C_27_to_66.pdb"
[atomselect top "chain C and protein and resid 80 to 95"] writepdb "C_80_to_95.pdb"
[atomselect top "chain C and protein and resid 99 to 140"] writepdb "C_99_to_140.pdb"
[atomselect top "chain C and protein and resid 157 to 176"] writepdb "C_157_to_176.pdb"
[atomselect top "chain C and protein and resid 187 to 245"] writepdb "C_187_to_245.pdb"
[atomselect top "chain C and protein and resid 261 to 443"] writepdb "C_261_to_443.pdb"
[atomselect top "chain C and protein and resid 449 to 454"] writepdb "C_449_to_454.pdb"
[atomselect top "chain C and protein and resid 460 to 471"] writepdb "C_460_to_471.pdb"
[atomselect top "chain C and protein and resid 487 to 498"] writepdb "C_487_to_498.pdb"
[atomselect top "chain C and protein and resid 503 to 620"] writepdb "C_503_to_620.pdb"
[atomselect top "chain C and protein and resid 641 to 672"] writepdb "C_641_to_672.pdb"
[atomselect top "chain C and protein and resid 687 to 811"] writepdb "C_687_to_811.pdb"
[atomselect top "chain C and protein and resid 815 to 828"] writepdb "C_815_to_828.pdb"
[atomselect top "chain C and protein and resid 853 to 1146"] writepdb "C_853_to_1146.pdb"
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
   pdb C_80_to_95.pdb
   residue 96 GLU C
   residue 97 LYS C
   residue 98 SER C
   pdb C_99_to_140.pdb
   residue 141 LEU C
   residue 142 GLY C
   residue 143 VAL C
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
   pdb C_157_to_176.pdb
   residue 177 MET C
   residue 178 ASP C
   residue 179 LEU C
   residue 180 GLU C
   residue 181 GLY C
   residue 182 LYS C
   residue 183 GLN C
   residue 184 GLY C
   residue 185 ASN C
   residue 186 PHE C
   pdb C_187_to_245.pdb
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
   pdb C_261_to_443.pdb
   residue 444 LYS C
   residue 445 VAL C
   residue 446 GLY C
   residue 447 GLY C
   residue 448 ASN C
   pdb C_449_to_454.pdb
   residue 455 LEU C
   residue 456 PHE C
   residue 457 ARG C
   residue 458 LYS C
   residue 459 SER C
   pdb C_460_to_471.pdb
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
   pdb C_487_to_498.pdb
   residue 499 PRO C
   residue 500 THR C
   residue 501 ASN C
   residue 502 GLY C
   pdb C_503_to_620.pdb
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
   pdb C_641_to_672.pdb
   residue 673 SER C
   residue 674 TYR C
   residue 675 GLN C
   residue 676 THR C
   residue 677 GLN C
   residue 678 THR C
   residue 679 ASN C
   residue 680 SER C
   residue 681 PRO C
   residue 682 GLY C
   residue 683 SER C
   residue 684 ALA C
   residue 685 SER C
   residue 686 SER C
   pdb C_687_to_811.pdb
   residue 812 PRO C
   residue 813 SER C
   residue 814 LYS C
   pdb C_815_to_828.pdb
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
   pdb C_853_to_1146.pdb
}
coordpdb C_27_to_66.pdb C
coordpdb C_80_to_95.pdb C
coordpdb C_99_to_140.pdb C
coordpdb C_157_to_176.pdb C
coordpdb C_187_to_245.pdb C
coordpdb C_261_to_443.pdb C
coordpdb C_449_to_454.pdb C
coordpdb C_460_to_471.pdb C
coordpdb C_487_to_498.pdb C
coordpdb C_503_to_620.pdb C
coordpdb C_641_to_672.pdb C
coordpdb C_687_to_811.pdb C
coordpdb C_815_to_828.pdb C
coordpdb C_853_to_1146.pdb C
coord C 67 N [cacoIn_nOut 66 C 0]
coord C 96 N [cacoIn_nOut 95 C 0]
coord C 141 N [cacoIn_nOut 140 C 0]
coord C 177 N [cacoIn_nOut 176 C 0]
coord C 246 N [cacoIn_nOut 245 C 0]
coord C 444 N [cacoIn_nOut 443 C 0]
coord C 455 N [cacoIn_nOut 454 C 0]
coord C 472 N [cacoIn_nOut 471 C 0]
coord C 499 N [cacoIn_nOut 498 C 0]
coord C 621 N [cacoIn_nOut 620 C 0]
coord C 673 N [cacoIn_nOut 672 C 0]
coord C 812 N [cacoIn_nOut 811 C 0]
coord C 829 N [cacoIn_nOut 828 C 0]
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
patch DISU B:336 B:361
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
patch DISU C:743 C:749
patch DISU C:1032 C:1043
patch DISU C:1082 C:1126
##### output of python3 parse_pdb_psfgen.py 6vsb.pdb above ###

guesscoord

regenerate angles dihedrals

writepsf "my_6vsb.psf"
writepdb "unrelaxed.pdb"

lappend LOCALFILES unrelaxed.pdb

mol delete top
mol new my_6vsb.psf
mol addfile unrelaxed.pdb
set molid [molinfo top get id]
set or [measure center [atomselect top "all"] weight mass]
set a [atomselect top all]
$a moveby [vecscale -1 $or]
if { $LOG_DCD != "0" } {
   mol new my_6vsb.psf
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
set temperature 2.5
set k 10.0
set r0 1.5
set bg [atomselect ${molid} "noh"]
foreach l $loops {
  set chain [lindex $l 0]
  set residueList [[atomselect ${molid} "chain $chain and resid [lindex $l 1] to [lindex $l 2] and name CA"] get residue]
  do_loop_mc ${residueList} ${chain} ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999] $logid
}
$a writepdb "my_6vsb_mcOut.pdb"

if { $LOG_DCD != "0" } {
   set loga [atomselect $logid all]
   animate write dcd $log_dcd_file waitfor all sel $loga $logid
}

# make a pdb file that fixes all heavy atoms in the original
# crystal structure -- all added atoms are set as unfixed
# for a minimization
mol delete top
mol new my_6vsb.psf
mol addfile unrelaxed.pdb
set a [atomselect top all]
$a set beta 0
foreach s $segs {
  [atomselect top "chain [lindex $s 0] and resid [lindex $s 1] to [lindex $s 2]"] set beta 1
}
[atomselect top "not noh"] set beta 0

$a writepdb "my_6vsb_fix.pdb"

# clean up
foreach f $LOCALFILES { 
  exec rm $f
}

quit

