# VMD script for generating solvated/neutralized system from the protein psf/pdb
#
# cameron f abrams (c) 2018
# cfa22@drexel.edu
# drexel university
# chemical and biological engineering

set seed 12345 
for {set i 0} {$i < $argc} {incr i} {
   set arg [lindex $argv i]
   if { $arg == "-seed" } {
      incr i
      set seed [lindex $argv $i]
   }
}
set LOCALFILES {}
set pad 35 

set inputname my_5jyn
set PSF ${inputname}.psf
set outputname1 ${inputname}_wb
set outputname2 ${inputname}_i

mol new $PSF
mol addfile ${inputname}_vac.coor

set a [atomselect top all]

set c [measure center $a]

$a moveby [vecscale $c -1.0]

set box { { ? ? ? } { ? ? ? } }
set basisvec { ? ? ? }
set origin { ? ? ? }

set minmax [measure minmax $a]
set diameterx [expr [lindex $minmax 1 0]-[lindex $minmax 0 0]]
set diametery [expr [lindex $minmax 1 1]-[lindex $minmax 0 1]]
set diameter $diameterx
if { [expr $diameter > $diametery] } {
   set diameter $diametery
}
set protein_area [expr ($diameter+4)*($diameter+4) * 3.141593 / 4.0]

foreach d {0 1} {
  lset box 0 $d [format "%.6f" [expr [lindex $minmax 0 $d] - $pad]]
  lset box 1 $d [format "%.6f" [expr [lindex $minmax 1 $d] + $pad]]
}
lset box 0 2 [expr [lindex $minmax 0 2] - 20]
lset box 1 2 [expr [lindex $minmax 1 2] + 20]

foreach d {0 1 2} {
  lset basisvec $d [expr [lindex $box 1 $d ] - [lindex $box 0 $d]] 
  lset origin $d [expr 0.5*([lindex $box 1 $d ] + [lindex $box 0 $d])] 
}

$a writepdb prot.pdb
lappend LOCALFILES prot.pdb

set mp [vecsum [$a get mass]]
set z [vecsum [$a get charge]]
set lx [format "%.6f" [expr [lindex [lindex $box 1] 0] - [lindex [lindex $box 0] 0]]]
set ly [format "%.6f" [expr [lindex [lindex $box 1] 1] - [lindex [lindex $box 0] 1]]]
set lz [format "%.6f" [expr [lindex [lindex $box 1] 2] - [lindex [lindex $box 0] 2]]]
set V [format "%.6f" [expr ($lx-2) * ($ly-2) * $lz]]
set A [format "%.6f" [expr ($lx-2) * ($ly-2)]]
set Vw [format "%.6f" [expr ($lx-2) * ($ly-2) * ($lz - 52. - 2)]]

set xmin [format "%.6f" [expr [lindex [lindex $box 0] 0]+1]]
set xmax [format "%.6f" [expr [lindex [lindex $box 1] 0]-1]]
set ymin [format "%.6f" [expr [lindex [lindex $box 0] 1]+1]]
set ymax [format "%.6f" [expr [lindex [lindex $box 1] 1]-1]]
set zmin [format "%.6f" [lindex [lindex $box 0] 2]]
set zmax [format "%.6f" [lindex [lindex $box 1] 2]]

set zmp [format "%.6f" [expr ($zmax + $zmin)/2.0]]
set zllower [format "%.6f" [expr $zmp - 29]]
set zlupper [format "%.6f" [expr $zmp + 29]]
set zwlower [format "%.6f" [expr $zllower - 2]]
set zwupper [format "%.6f" [expr $zlupper + 2]]

puts "$mp $z $lx $ly $lz $V $A $diameter $protein_area"

# SAPL for DMPC is 60.6 A^2 (Nagle, Biophys J, 2005; 10.1529/biophysj.104.056606)
# but using 64 per charmm-gui
set SAPLFAC 1.33
set SAPL [expr 64.0*$SAPLFAC]
# assume protein complex's XY-projection is circular 
set AvailA [expr $A - $protein_area]
set nLipid [expr int($AvailA/$SAPL)]
#set nLipid 5

set MWw 18.0
set densgcc 1.0
set nw [expr int( 1.0 / $MWw * (0.6022*$densgcc*$Vw) )]
#set nw 200

set nna 0
set ncl 0
if { $z > 0 } {
   set ncl [expr round($z)]
} elseif { $z < 0 } {
   set nna [expr round(-($z))]
}

puts "system will have $nw waters, [expr 2 * $nLipid] lipids, $nna sodiums and $ncl chlorides."

set fp [open "pm-tmp.in" "w"]
puts $fp "
output my_5jyn_packed.pdb
filetype pdb
seed $seed
tolerance 2.0
structure prot.pdb
  number 1
  center
  fixed 0. 0. 0. 0. 0. 0.
end structure
# z(-) leaflet
structure dmpc.pdb
  number $nLipid
  inside box $xmin $ymin $zllower $xmax $ymax $zmp
  atoms 1
    below plane 0. 0. 1. [expr $zllower + 2]
  end atoms
  atoms 114
    over plane 0. 0. 1. [expr $zmp - 2]
  end atoms
end structure
# z(+) leaflet
structure dmpc.pdb
  number $nLipid
  inside box $xmin $ymin $zmp $xmax $ymax $zlupper
  atoms 1
    over plane 0. 0. 1. [expr $zlupper - 2]
  end atoms
  atoms 114
    below plane 0. 0. 1. [expr $zmp + 2]
  end atoms
end structure
structure TIP3.pdb
  resnumbers 0
  number [expr $nw / 2]
  inside box $xmin $ymin $zmin $xmax $ymax $zwlower
end structure
structure TIP3.pdb
  resnumbers 0
  number [expr $nw / 2]
  inside box $xmin $ymin $zwupper $xmax $ymax $zmax
end structure
"

if { $nna > 0 } {
   set nna_lower [expr $nna / 2]
   set nna_upper [expr $nna - $nna_lower]
   puts $fp "
structure SOD.pdb
   resnumbers 0
   number $nna_lower
   inside box $xmin $ymin $zmin $xmax $ymax $zwlower
end structure
"
   if { $nna_upper > 0 } {
      puts $fp "
structure SOD.pdb
   resnumbers 0
   number $nna_upper
   inside box $xmin $ymin $zwupper $xmax $ymax $zmax
end structure
"
  }
}

if { $ncl > 0 } {
   set ncl_lower [expr $ncl / 2]
   set ncl_upper [expr $ncl - $ncl_lower]
   puts $fp "
structure CL.pdb
   resnumbers 0
   number $ncl_lower
   inside box $xmin $ymin $zmin $xmax $ymax $zwlower
end structure
"
   if { $ncl_upper > 0 } {
      puts $fp "
structure CL.pdb
   resnumbers 0
   number $ncl_upper
   inside box $xmin $ymin $zwupper $xmax $ymax $zmax
end structure
" 
   }
}
close $fp

set water_file [ open TIP3.pdb w ]
puts $water_file "
HETATM    1  H1  TIP3    1       9.626   6.787  12.673  0.00  0.00      W
HETATM    2  H2  TIP3    1       9.626   8.420  12.673  0.00  0.00      W
HETATM    3  OH2 TIP3    1      10.203   7.604  12.673  0.00  0.00      W"
close $water_file

set sod_file [ open SOD.pdb w ]
puts $sod_file  "HETATM    1 SOD  SOD    2        0.000   0.000   0.000  0.00  0.00      I" 
close $sod_file

set cl_file [ open CL.pdb w ] 
puts $cl_file  "HETATM    1 CLA  CLA    2        0.000   0.000   0.000  0.00  0.00      I" 
close $cl_file 

set lipid_file [ open dmpc.pdb w ]
# the following is a pdbfile generated by running psfgen/guesscoord on 
# positions for N, C13, and C14 of DMPC, followed by rotation of one
# C-C bond to bring the two tails into a parallel arrangement
puts $lipid_file "
ATOM      1  N   DMPCA  62     -30.964 -26.029  17.702  1.00  0.00      A    N
ATOM      2  C13 DMPCA  62     -30.069 -26.317  18.925  1.00  0.00      A    C
ATOM      3 H13A DMPCA  62     -30.630 -26.797  19.715  0.00  0.00      A    H
ATOM      4 H13B DMPCA  62     -29.264 -26.986  18.633  0.00  0.00      A    H
ATOM      5 H13C DMPCA  62     -29.621 -25.410  19.307  0.00  0.00      A    H
ATOM      6  C14 DMPCA  62     -31.598 -27.410  17.221  1.00  0.00      A    C
ATOM      7 H14A DMPCA  62     -30.830 -28.072  16.807  0.00  0.00      A    H
ATOM      8 H14B DMPCA  62     -32.140 -27.930  17.989  0.00  0.00      A    H
ATOM      9 H14C DMPCA  62     -32.249 -27.247  16.357  0.00  0.00      A    H
ATOM     10  C15 DMPCA  62     -32.130 -25.194  18.232  1.00  0.00      A    C
ATOM     11 H15A DMPCA  62     -32.646 -25.714  19.029  0.00  0.00      A    H
ATOM     12 H15B DMPCA  62     -31.794 -24.235  18.601  0.00  0.00      A    H
ATOM     13 H15C DMPCA  62     -32.847 -25.019  17.434  0.00  0.00      A    H
ATOM     14  C12 DMPCA  62     -30.235 -25.406  16.530  1.00  0.00      A    C
ATOM     15 H12A DMPCA  62     -29.923 -24.395  16.741  0.00  0.00      A    H
ATOM     16 H12B DMPCA  62     -30.936 -25.335  15.689  0.00  0.00      A    H
ATOM     17  C11 DMPCA  62     -29.035 -26.188  15.961  0.00  0.00      A    C
ATOM     18 H11A DMPCA  62     -28.887 -27.180  16.445  0.00  0.00      A    H
ATOM     19 H11B DMPCA  62     -28.103 -25.606  16.137  0.00  0.00      A    H
ATOM     20  P   DMPCA  62     -30.405 -27.265  14.048  0.00  0.00      A    P
ATOM     21  O13 DMPCA  62     -30.326 -28.551  14.774  0.00  0.00      A    O
ATOM     22  O14 DMPCA  62     -31.660 -26.501  14.240  0.00  0.00      A    O
ATOM     23  O11 DMPCA  62     -30.101 -27.468  12.513  0.00  0.00      A    O
ATOM     24  O12 DMPCA  62     -29.222 -26.337  14.558  0.00  0.00      A    O
ATOM     25  C1  DMPCA  62     -30.968 -28.322  11.759  0.00  0.00      A    C
ATOM     26  HA  DMPCA  62     -32.043 -28.074  11.940  0.00  0.00      A    H
ATOM     27  HB  DMPCA  62     -30.792 -29.387  12.039  0.00  0.00      A    H
ATOM     28  C2  DMPCA  62     -30.715 -28.163  10.238  0.00  0.00      A    C
ATOM     29  HS  DMPCA  62     -29.645 -28.431  10.061  0.00  0.00      A    H
ATOM     30  O21 DMPCA  62     -31.002 -26.797   9.877  0.00  0.00      A    O
ATOM     31  C21 DMPCA  62     -30.270 -26.314   8.893  0.00  0.00      A    C
ATOM     32  O22 DMPCA  62     -29.405 -26.907   8.274  0.00  0.00      A    O
ATOM     33  C22 DMPCA  62     -30.698 -24.871   8.625  0.00  0.00      A    C
ATOM     34  H2R DMPCA  62     -30.441 -24.270   9.521  0.00  0.00      A    H
ATOM     35  H2S DMPCA  62     -31.799 -24.862   8.493  0.00  0.00      A    H
ATOM     36  C3  DMPCA  62     -31.620 -29.137   9.426  0.00  0.00      A    C
ATOM     37  HX  DMPCA  62     -32.671 -28.909   9.714  0.00  0.00      A    H
ATOM     38  HY  DMPCA  62     -31.419 -30.185   9.746  0.00  0.00      A    H
ATOM     39  O31 DMPCA  62     -31.497 -28.953   7.999  0.00  0.00      A    O
ATOM     40  C31 DMPCA  62     -30.516 -29.645   7.423  0.00  0.00      A    C
ATOM     41  O32 DMPCA  62     -29.840 -30.501   7.963  0.00  0.00      A    O
ATOM     42  C32 DMPCA  62     -30.402 -29.198   5.966  0.00  0.00      A    C
ATOM     43  H2X DMPCA  62     -30.014 -28.157   5.971  0.00  0.00      A    H
ATOM     44  H2Y DMPCA  62     -31.420 -29.194   5.528  0.00  0.00      A    H
ATOM     45  C23 DMPCA  62     -30.009 -24.287   7.371  0.00  0.00      A    C
ATOM     46  H3R DMPCA  62     -28.908 -24.245   7.537  0.00  0.00      A    H
ATOM     47  H3S DMPCA  62     -30.361 -23.243   7.201  0.00  0.00      A    H
ATOM     48  C24 DMPCA  62     -30.282 -25.117   6.111  0.00  0.00      A    C
ATOM     49  H4R DMPCA  62     -31.376 -25.137   5.910  0.00  0.00      A    H
ATOM     50  H4S DMPCA  62     -29.949 -26.166   6.286  0.00  0.00      A    H
ATOM     51  C25 DMPCA  62     -29.551 -24.563   4.880  0.00  0.00      A    C
ATOM     52  H5R DMPCA  62     -28.458 -24.526   5.096  0.00  0.00      A    H
ATOM     53  H5S DMPCA  62     -29.891 -23.520   4.690  0.00  0.00      A    H
ATOM     54  C26 DMPCA  62     -29.781 -25.408   3.620  0.00  0.00      A    C
ATOM     55  H6R DMPCA  62     -30.867 -25.426   3.383  0.00  0.00      A    H
ATOM     56  H6S DMPCA  62     -29.460 -26.455   3.822  0.00  0.00      A    H
ATOM     57  C27 DMPCA  62     -29.005 -24.873   2.409  0.00  0.00      A    C
ATOM     58  H7R DMPCA  62     -27.921 -24.840   2.660  0.00  0.00      A    H
ATOM     59  H7S DMPCA  62     -29.333 -23.832   2.194  0.00  0.00      A    H
ATOM     60  C28 DMPCA  62     -29.199 -25.734   1.154  0.00  0.00      A    C
ATOM     61  H8R DMPCA  62     -30.278 -25.752   0.884  0.00  0.00      A    H
ATOM     62  H8S DMPCA  62     -28.889 -26.779   1.380  0.00  0.00      A    H
ATOM     63  C29 DMPCA  62     -28.385 -25.220  -0.041  0.00  0.00      A    C
ATOM     64  H9R DMPCA  62     -27.308 -25.193   0.238  0.00  0.00      A    H
ATOM     65  H9S DMPCA  62     -28.698 -24.179  -0.276  0.00  0.00      A    H
ATOM     66 C210 DMPCA  62     -28.555 -26.093  -1.292  0.00  0.00      A    C
ATOM     67 H10R DMPCA  62     -29.628 -26.107  -1.585  0.00  0.00      A    H
ATOM     68 H10S DMPCA  62     -28.256 -27.138  -1.049  0.00  0.00      A    H
ATOM     69 C211 DMPCA  62     -27.714 -25.596  -2.474  0.00  0.00      A    C
ATOM     70 H11R DMPCA  62     -26.642 -25.576  -2.175  0.00  0.00      A    H
ATOM     71 H11S DMPCA  62     -28.015 -24.552  -2.724  0.00  0.00      A    H
ATOM     72 C212 DMPCA  62     -27.869 -26.477  -3.722  0.00  0.00      A    C
ATOM     73 H12R DMPCA  62     -28.939 -26.490  -4.027  0.00  0.00      A    H
ATOM     74 H12S DMPCA  62     -27.577 -27.520  -3.468  0.00  0.00      A    H
ATOM     75 C213 DMPCA  62     -27.016 -25.989  -4.899  0.00  0.00      A    C
ATOM     76 H13R DMPCA  62     -25.947 -25.970  -4.584  0.00  0.00      A    H
ATOM     77 H13S DMPCA  62     -27.310 -24.944  -5.147  0.00  0.00      A    H
ATOM     78 C214 DMPCA  62     -27.158 -26.862  -6.149  0.00  0.00      A    C
ATOM     79 H14R DMPCA  62     -26.848 -27.907  -5.931  0.00  0.00      A    H
ATOM     80 H14S DMPCA  62     -28.213 -26.876  -6.496  0.00  0.00      A    H
ATOM     81 H14T DMPCA  62     -26.522 -26.473  -6.973  0.00  0.00      A    H
ATOM     82  C33 DMPCA  62     -29.469 -30.094   5.120  0.00  0.00      A    C
ATOM     83  H3X DMPCA  62     -29.842 -31.143   5.136  0.00  0.00      A    H
ATOM     84  H3Y DMPCA  62     -28.445 -30.093   5.558  0.00  0.00      A    H
ATOM     85  C34 DMPCA  62     -29.397 -29.609   3.667  0.00  0.00      A    C
ATOM     86  H4X DMPCA  62     -28.980 -28.576   3.649  0.00  0.00      A    H
ATOM     87  H4Y DMPCA  62     -30.429 -29.563   3.251  0.00  0.00      A    H
ATOM     88  C35 DMPCA  62     -28.544 -30.511   2.766  0.00  0.00      A    C
ATOM     89  H5X DMPCA  62     -28.953 -31.546   2.797  0.00  0.00      A    H
ATOM     90  H5Y DMPCA  62     -27.502 -30.546   3.156  0.00  0.00      A    H
ATOM     91  C36 DMPCA  62     -28.529 -30.020   1.312  0.00  0.00      A    C
ATOM     92  H6X DMPCA  62     -28.084 -29.000   1.275  0.00  0.00      A    H
ATOM     93  H6Y DMPCA  62     -29.578 -29.944   0.948  0.00  0.00      A    H
ATOM     94  C37 DMPCA  62     -27.748 -30.945   0.369  0.00  0.00      A    C
ATOM     95  H7X DMPCA  62     -28.181 -31.969   0.426  0.00  0.00      A    H
ATOM     96  H7Y DMPCA  62     -26.688 -31.006   0.705  0.00  0.00      A    H
ATOM     97  C38 DMPCA  62     -27.795 -30.462  -1.086  0.00  0.00      A    C
ATOM     98  H8X DMPCA  62     -27.339 -29.448  -1.150  0.00  0.00      A    H
ATOM     99  H8Y DMPCA  62     -28.859 -30.373  -1.402  0.00  0.00      A    H
ATOM    100  C39 DMPCA  62     -27.069 -31.404  -2.055  0.00  0.00      A    C
ATOM    101  H9X DMPCA  62     -27.515 -32.422  -1.975  0.00  0.00      A    H
ATOM    102  H9Y DMPCA  62     -25.998 -31.480  -1.761  0.00  0.00      A    H
ATOM    103 C310 DMPCA  62     -27.165 -30.928  -3.510  0.00  0.00      A    C
ATOM    104 H10X DMPCA  62     -26.709 -29.915  -3.593  0.00  0.00      A    H
ATOM    105 H10Y DMPCA  62     -28.238 -30.839  -3.792  0.00  0.00      A    H
ATOM    106 C311 DMPCA  62     -26.468 -31.876  -4.494  0.00  0.00      A    C
ATOM    107 H11X DMPCA  62     -26.918 -32.890  -4.404  0.00  0.00      A    H
ATOM    108 H11Y DMPCA  62     -25.391 -31.958  -4.223  0.00  0.00      A    H
ATOM    109 C312 DMPCA  62     -26.588 -31.401  -5.948  0.00  0.00      A    C
ATOM    110 H12X DMPCA  62     -26.137 -30.387  -6.040  0.00  0.00      A    H
ATOM    111 H12Y DMPCA  62     -27.665 -31.316  -6.216  0.00  0.00      A    H
ATOM    112 C313 DMPCA  62     -25.899 -32.348  -6.938  0.00  0.00      A    C
ATOM    113 H13X DMPCA  62     -26.349 -33.362  -6.836  0.00  0.00      A    H
ATOM    114 H13Y DMPCA  62     -24.822 -32.431  -6.664  0.00  0.00      A    H
ATOM    115 C314 DMPCA  62     -26.014 -31.890  -8.395  0.00  0.00      A    C
ATOM    116 H14X DMPCA  62     -27.081 -31.823  -8.697  0.00  0.00      A    H
ATOM    117 H14Y DMPCA  62     -25.550 -30.888  -8.525  0.00  0.00      A    H
ATOM    118 H14Z DMPCA  62     -25.501 -32.606  -9.071  0.00  0.00      A    H
"
close $lipid_file

# generate an input file for the first solvated MD simulation
# namd config file
set fp [open "cell.inp" "w"]
puts $fp "cellbasisvector1 [lindex $basisvec 0] 0 0"
puts $fp "cellbasisvector2 0 [lindex $basisvec 1] 0"
puts $fp "cellbasisvector3 0 0 [lindex $basisvec 2]"
puts $fp "cellorigin $origin"
close $fp
puts "Generated cell.inp."

exit
