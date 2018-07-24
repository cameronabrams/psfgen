# VMD script for generating solvated/neutralized system from the protein psf/pdb
#
# cameron f abrams (c) 2018
# cfa22@drexel.edu
# drexel university
# chemical and biological engineering

set inputname my_5jyn
set seed 12345
# pad in Angstoms along +/-X,Y from the protein to define the lateral box size 
set xypad 35.0
# pad in Angstroms along +/-Z from the protein to define the box height
set zpad 20.0
# initial depth of lipid leaflet from N to C-tail; referenced to a all-trans DMPC
set leaflet_depth 30.0
# factor for scaling SAPL; less than one since lipid template is all-trans
set SAPLFAC 0.90
# density of water layer in g/cc
set watergcc 1.0
for { set i 0 } { $i < $argc } { incr  i} {
   set arg [lindex $argv $i]
   if { $arg == "-seed" } {
      incr i
      set seed [lindex $argv $i]
   }
   if { $arg == "-xypad" } {
      incr i
      set xypad [lindex $argv $i]
   }
   if { $arg == "-zpad" } {
      incr i
      set zpad [lindex $argv $i]
   }
   if { $arg == "-leaflet_depth" } {
      incr i
      set leaflet_depth [lindex $argv $i]
   }
   if { $arg == "-saplfac" } {
      incr i
      set SAPLFAC [lindex $argv $i]
   }
   if { $arg == "-watergcc" } {
      incr i
      set watergcc [lindex $argv $i]
   }
}
set LOCALFILES {}

set PSF ${inputname}.psf
set outputname1 ${inputname}_wb
set outputname2 ${inputname}_i

mol new $PSF
mol addfile ${inputname}_vac.coor
set a [atomselect top all]
set mp [vecsum [$a get mass]]
set z [vecsum [$a get charge]]
set c [measure center $a]
$a moveby [vecscale $c -1.0]
$a writepdb prot.pdb; lappend LOCALFILES prot.pdb
set minmax [measure minmax $a]
# use a measurement of the moment of inertia to determine
# the cross-sectional area of the protein
set ev1 [lindex [lindex [measure inertia $a eigenvals] end] 0]
set radius [expr sqrt( $ev1 / $mp )]
set protein_area [expr ($radius)*($radius) * 3.141593 ]
# done with molecule data
mol delete top

set box { { ? ? ? } { ? ? ? } }
set basisvec { ? ? ? }
set origin { ? ? ? }

# compute box lower left and upper right corners 
# using pad values provided by user or defaults
foreach d {0 1} {
  lset box 0 $d [format "%.6f" [expr [lindex $minmax 0 $d] - $xypad]]
  lset box 1 $d [format "%.6f" [expr [lindex $minmax 1 $d] + $xypad]]
}
lset box 0 2 [expr [lindex $minmax 0 2] - $zpad]
lset box 1 2 [expr [lindex $minmax 1 2] + $zpad]

foreach d {0 1 2} {
  lset basisvec $d [expr [lindex $box 1 $d ] - [lindex $box 0 $d]] 
  lset origin $d [expr 0.5*([lindex $box 1 $d ] + [lindex $box 0 $d])] 
}

# compute box dimensions
set lx [format "%.6f" [expr [lindex [lindex $box 1] 0] - [lindex [lindex $box 0] 0]]]
set ly [format "%.6f" [expr [lindex [lindex $box 1] 1] - [lindex [lindex $box 0] 1]]]
set lz [format "%.6f" [expr [lindex [lindex $box 1] 2] - [lindex [lindex $box 0] 2]]]

# save corner coordinates of box into new variables; pad by 1 A for PBC (only in Z)
set xmin [format "%.6f" [expr [lindex [lindex $box 0] 0]]]
set xmax [format "%.6f" [expr [lindex [lindex $box 1] 0]]]
set ymin [format "%.6f" [expr [lindex [lindex $box 0] 1]]]
set ymax [format "%.6f" [expr [lindex [lindex $box 1] 1]]]
set zmin [format "%.6f" [expr [lindex [lindex $box 0] 2]+1]]
set zmax [format "%.6f" [expr [lindex [lindex $box 1] 2]-1]]

# compute cross-sectional area A and volume V
set A [format "%.6f" [expr $lx * $ly]]
set V [format "%.6f" [expr $A * $lz]]

# midplane z-position
set zmp [format "%.6f" [expr ($zmax + $zmin)/2.0]]
# z-coordinate of upper boundary of lower leaflet
set zLLhi $zmp
# z-coordinate of lower boundary of upper leaflet
set zULlo $zmp
# z-coordinate of lower boundary of lower leaflet
set zLLlo [format "%.6f" [expr $zmp - $leaflet_depth]] 
# z-coordinate of upper boundary of upper leaflet
set zULhi [format "%.6f" [expr $zmp + $leaflet_depth]]
# z-coordinate of upper boundary of lower water layer
set zLWhi [format "%.6f" [expr $zLLlo - 1]]
# z-coordinate of lower boundary of lower water layer
set zLWlo $zmin
# z-coordinate of lower boundary or upper water layer
set zUWlo [format "%.6f" [expr $zULhi + 1]]
# z-coordinate of upper boundary of upper water layer
set zUWhi $zmax

# water volume
set Vw [format "%.6f" [expr $A * (($zLWhi - $zLWlo) + ($zUWhi - $zUWlo))]]

# SAPL for DMPC is 60.6 A^2 (Nagle, Biophys J, 2005; 10.1529/biophysj.104.056606)
set SAPL [expr 60.6*$SAPLFAC]
# assume protein complex's XY-projection is circular 
# compute available cross-sectional area for packing lipids
set AvailA [expr $A - $protein_area]
# compute the number of lipids per leaflet
set nLipid [expr int($AvailA/$SAPL)]
# compute the TOTAL number of waters in BOTH layers
set nw [expr int( 1.0 / 18.0 * (0.6022*$watergcc*$Vw) )]

# compute the number of solvent ions necessary to neutralize only
# this is the total in both layers
set nna 0
set ncl 0
if { $z > 0 } {
   set ncl [expr round($z)]
} elseif { $z < 0 } {
   set nna [expr round(-($z))]
}

puts "--------------- [format "%12.6f"  $zmax] -----------------"
puts "         waters: [expr int($nw/2)]"
puts "             cl: [expr int($ncl/2)]"
puts "             na: [expr int($nna/2)]"
puts "--------------- [format "%12.6f" $zUWlo] -----------------"
puts "--------------- [format "%12.6f" $zULhi] -----------------"
puts "                          O"
puts "                         /|"
puts "                        / | x $nLipid; gap size [expr $zULhi - $zULlo] A"  
puts "                        | |"
puts "                        | |"
puts "--------------- [format "%12.6f" $zULlo] -----------------"
puts "----------------[format "%12.6f"   $zmp] -----------------"
puts "--------------- [format "%12.6f" $zLLhi] -----------------"
puts "                          | |"
puts "                          | | x $nLipid; gap size [expr $zLLhi - $zLLlo] A"
puts "                          | /"
puts "                          |/"
puts "                          O"
puts "--------------- [format "%12.6f" $zLLlo] -----------------"
puts "--------------- [format "%12.6f" $zLWhi] -----------------"
puts "         waters: [expr int($nw/2)]"
puts "             cl: [expr int($ncl/2)]"
puts "             na: [expr int($nna/2)]"
puts "--------------- [format "%12.6f"  $zmin] -----------------"

# generate packmol input files
# stage 1 packs the protein at the middle and then packs just the lower leaflet
set fp [open "pm-stage1.in" "w"]
puts $fp "
output my_5jyn_packed_stage1.pdb
filetype pdb
seed $seed
tolerance 2.0
structure prot.pdb
  number 1
  fixed 0. 0. 0. 0. 0. 0.
end structure
# z(-) leaflet
structure dmpc.pdb
  number $nLipid
  inside box $xmin $ymin $zLLlo $xmax $ymax $zLLhi
  constrain_rotation x 180.0 1.0
  constrain_rotation y 180.0 1.0
end structure
"
close $fp
# stage 2 packs only the upper leaflet
set fp [open "pm-stage2.in" "w"]
puts $fp "
output my_5jyn_packed_stage2.pdb
filetype pdb
seed $seed
tolerance 2.0
structure my_5jyn_packed_stage1.pdb
  number 1
  fixed 0. 0. 0. 0. 0. 0.
end structure
# z(+) leaflet
structure dmpc.pdb
  number $nLipid
  inside box $xmin $ymin $zULlo $xmax $ymax $zULhi
  constrain_rotation x 0.0 1.0
  constrain_rotation y 0.0 1.0
end structure
"
close $fp
# stage 3 packs the waters and solvent ions
set fp [open "pm-stage3.in" "w"]
puts $fp "
output my_5jyn_packed.pdb
filetype pdb
seed $seed
tolerance 2.0
structure my_5jyn_packed_stage2.pdb
  number 1
  fixed 0. 0. 0. 0. 0. 0.
end structure
structure TIP3.pdb
  resnumbers 0
  number [expr $nw / 2]
  inside box $xmin $ymin $zmin $xmax $ymax $zLWhi
end structure
structure TIP3.pdb
  resnumbers 0
  number [expr $nw / 2]
  inside box $xmin $ymin $zUWlo $xmax $ymax $zmax
end structure
"

if { $nna > 0 } {
   set nna_lower [expr $nna / 2]
   set nna_upper [expr $nna - $nna_lower]
   puts $fp "
structure SOD.pdb
   resnumbers 0
   number $nna_lower
   inside box $xmin $ymin $zmin $xmax $ymax $zLWhi
end structure
"
   if { $nna_upper > 0 } {
      puts $fp "
structure SOD.pdb
   resnumbers 0
   number $nna_upper
   inside box $xmin $ymin $zUWlo $xmax $ymax $zmax
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
   inside box $xmin $ymin $zmin $xmax $ymax $zLWhi
end structure
"
   if { $ncl_upper > 0 } {
      puts $fp "
structure CL.pdb
   resnumbers 0
   number $ncl_upper
   inside box $xmin $ymin $zUWlo $xmax $ymax $zmax
end structure
" 
   }
}
close $fp

# here we generate local copies of the template PDB files for TIP3P water, solvent ions,
# and a perfectly all-trans, verically oriented (head group high) DMPC molecule
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
# C-C bond to bring the two tails into a parallel arrangement, and 
# euler rotations to make the molecule as "vertical" as possible
puts $lipid_file "
ATOM      1  N   DMPCA  62       0.912  -0.508  13.425  1.00  0.00      A    N
ATOM      2  C13 DMPCA  62       2.024  -1.012  14.369  1.00  0.00      A    C
ATOM      3 H13A DMPCA  62       1.615  -1.617  15.166  0.00  0.00      A    H
ATOM      4 H13B DMPCA  62       2.723  -1.625  13.807  0.00  0.00      A    H
ATOM      5 H13C DMPCA  62       2.579  -0.188  14.797  0.00  0.00      A    H
ATOM      6  C14 DMPCA  62       0.135  -1.779  12.860  1.00  0.00      A    C
ATOM      7 H14A DMPCA  62       0.773  -2.365  12.189  0.00  0.00      A    H
ATOM      8 H14B DMPCA  62      -0.262  -2.422  13.625  0.00  0.00      A    H
ATOM      9 H14C DMPCA  62      -0.670  -1.463  12.191  0.00  0.00      A    H
ATOM     10  C15 DMPCA  62      -0.084   0.230  14.319  1.00  0.00      A    C
ATOM     11 H15A DMPCA  62      -0.449  -0.417  15.107  0.00  0.00      A    H
ATOM     12 H15B DMPCA  62       0.361   1.107  14.767  0.00  0.00      A    H
ATOM     13 H15C DMPCA  62      -0.940   0.547  13.729  0.00  0.00      A    H
ATOM     14  C12 DMPCA  62       1.414   0.304  12.250  1.00  0.00      A    C
ATOM     15 H12A DMPCA  62       1.806   1.260  12.559  0.00  0.00      A    H
ATOM     16 H12B DMPCA  62       0.561   0.526  11.598  0.00  0.00      A    H
ATOM     17  C11 DMPCA  62       2.438  -0.375  11.320  0.00  0.00      A    C
ATOM     18 H11A DMPCA  62       2.639  -1.437  11.588  0.00  0.00      A    H
ATOM     19 H11B DMPCA  62       3.411   0.160  11.395  0.00  0.00      A    H
ATOM     20  P   DMPCA  62       0.663  -1.092   9.579  0.00  0.00      A    P
ATOM     21  O13 DMPCA  62       0.832  -2.486  10.045  0.00  0.00      A    O
ATOM     22  O14 DMPCA  62      -0.493  -0.364  10.155  0.00  0.00      A    O
ATOM     23  O11 DMPCA  62       0.639  -1.026   8.004  0.00  0.00      A    O
ATOM     24  O12 DMPCA  62       1.964  -0.275   9.982  0.00  0.00      A    O
ATOM     25  C1  DMPCA  62      -0.399  -1.729   7.313  0.00  0.00      A    C
ATOM     26  HA  DMPCA  62      -1.403  -1.509   7.752  0.00  0.00      A    H
ATOM     27  HB  DMPCA  62      -0.216  -2.828   7.366  0.00  0.00      A    H
ATOM     28  C2  DMPCA  62      -0.454  -1.309   5.822  0.00  0.00      A    C
ATOM     29  HS  DMPCA  62       0.545  -1.551   5.383  0.00  0.00      A    H
ATOM     30  O21 DMPCA  62      -0.749   0.101   5.765  0.00  0.00      A    O
ATOM     31  C21 DMPCA  62      -0.212   0.743   4.746  0.00  0.00      A    C
ATOM     32  O22 DMPCA  62       0.482   0.261   3.870  0.00  0.00      A    O
ATOM     33  C22 DMPCA  62      -0.623   2.213   4.821  0.00  0.00      A    C
ATOM     34  H2R DMPCA  62      -0.164   2.647   5.732  0.00  0.00      A    H
ATOM     35  H2S DMPCA  62      -1.727   2.254   4.923  0.00  0.00      A    H
ATOM     36  C3  DMPCA  62      -1.547  -2.120   5.062  0.00  0.00      A    C
ATOM     37  HX  DMPCA  62      -2.506  -1.938   5.597  0.00  0.00      A    H
ATOM     38  HY  DMPCA  62      -1.330  -3.210   5.152  0.00  0.00      A    H
ATOM     39  O31 DMPCA  62      -1.708  -1.691   3.694  0.00  0.00      A    O
ATOM     40  C31 DMPCA  62      -0.896  -2.280   2.817  0.00  0.00      A    C
ATOM     41  O32 DMPCA  62      -0.162  -3.221   3.052  0.00  0.00      A    O
ATOM     42  C32 DMPCA  62      -1.061  -1.585   1.466  0.00  0.00      A    C
ATOM     43  H2X DMPCA  62      -0.636  -0.564   1.565  0.00  0.00      A    H
ATOM     44  H2Y DMPCA  62      -2.146  -1.498   1.255  0.00  0.00      A    H
ATOM     45  C23 DMPCA  62      -0.179   3.002   3.568  0.00  0.00      A    C
ATOM     46  H3R DMPCA  62       0.933   3.007   3.506  0.00  0.00      A    H
ATOM     47  H3S DMPCA  62      -0.513   4.062   3.654  0.00  0.00      A    H
ATOM     48  C24 DMPCA  62      -0.738   2.406   2.270  0.00  0.00      A    C
ATOM     49  H4R DMPCA  62      -1.850   2.430   2.301  0.00  0.00      A    H
ATOM     50  H4S DMPCA  62      -0.422   1.341   2.192  0.00  0.00      A    H
ATOM     51  C25 DMPCA  62      -0.250   3.161   1.026  0.00  0.00      A    C
ATOM     52  H5R DMPCA  62       0.865   3.152   1.013  0.00  0.00      A    H
ATOM     53  H5S DMPCA  62      -0.576   4.224   1.090  0.00  0.00      A    H
ATOM     54  C26 DMPCA  62      -0.767   2.551  -0.283  0.00  0.00      A    C
ATOM     55  H6R DMPCA  62      -1.878   2.582  -0.288  0.00  0.00      A    H
ATOM     56  H6S DMPCA  62      -0.458   1.482  -0.332  0.00  0.00      A    H
ATOM     57  C27 DMPCA  62      -0.231   3.283  -1.520  0.00  0.00      A    C
ATOM     58  H7R DMPCA  62       0.881   3.264  -1.498  0.00  0.00      A    H
ATOM     59  H7S DMPCA  62      -0.551   4.348  -1.483  0.00  0.00      A    H
ATOM     60  C28 DMPCA  62      -0.714   2.655  -2.835  0.00  0.00      A    C
ATOM     61  H8R DMPCA  62      -1.825   2.693  -2.873  0.00  0.00      A    H
ATOM     62  H8S DMPCA  62      -0.410   1.585  -2.858  0.00  0.00      A    H
ATOM     63  C29 DMPCA  62      -0.138   3.364  -4.069  0.00  0.00      A    C
ATOM     64  H9R DMPCA  62       0.973   3.334  -4.019  0.00  0.00      A    H
ATOM     65  H9S DMPCA  62      -0.447   4.433  -4.054  0.00  0.00      A    H
ATOM     66 C210 DMPCA  62      -0.597   2.724  -5.386  0.00  0.00      A    C
ATOM     67 H10R DMPCA  62      -1.706   2.769  -5.447  0.00  0.00      A    H
ATOM     68 H10S DMPCA  62      -0.300   1.650  -5.390  0.00  0.00      A    H
ATOM     69 C211 DMPCA  62       0.007   3.414  -6.616  0.00  0.00      A    C
ATOM     70 H11R DMPCA  62       1.117   3.373  -6.547  0.00  0.00      A    H
ATOM     71 H11S DMPCA  62      -0.293   4.487  -6.617  0.00  0.00      A    H
ATOM     72 C212 DMPCA  62      -0.436   2.764  -7.935  0.00  0.00      A    C
ATOM     73 H12R DMPCA  62      -1.546   2.813  -8.009  0.00  0.00      A    H
ATOM     74 H12S DMPCA  62      -0.144   1.691  -7.927  0.00  0.00      A    H
ATOM     75 C213 DMPCA  62       0.180   3.444  -9.163  0.00  0.00      A    C
ATOM     76 H13R DMPCA  62       1.290   3.401  -9.078  0.00  0.00      A    H
ATOM     77 H13S DMPCA  62      -0.113   4.518  -9.165  0.00  0.00      A    H
ATOM     78 C214 DMPCA  62      -0.251   2.803 -10.486  0.00  0.00      A    C
ATOM     79 H14R DMPCA  62       0.051   1.734 -10.517  0.00  0.00      A    H
ATOM     80 H14S DMPCA  62      -1.354   2.858 -10.603  0.00  0.00      A    H
ATOM     81 H14T DMPCA  62       0.220   3.325 -11.346  0.00  0.00      A    H
ATOM     82  C33 DMPCA  62      -0.360  -2.327   0.305  0.00  0.00      A    C
ATOM     83  H3X DMPCA  62      -0.766  -3.360   0.221  0.00  0.00      A    H
ATOM     84  H3Y DMPCA  62       0.731  -2.410   0.515  0.00  0.00      A    H
ATOM     85  C34 DMPCA  62      -0.564  -1.596  -1.028  0.00  0.00      A    C
ATOM     86  H4X DMPCA  62      -0.115  -0.579  -0.957  0.00  0.00      A    H
ATOM     87  H4Y DMPCA  62      -1.656  -1.472  -1.207  0.00  0.00      A    H
ATOM     88  C35 DMPCA  62       0.048  -2.334  -2.226  0.00  0.00      A    C
ATOM     89  H5X DMPCA  62      -0.390  -3.355  -2.286  0.00  0.00      A    H
ATOM     90  H5Y DMPCA  62       1.145  -2.444  -2.072  0.00  0.00      A    H
ATOM     91  C36 DMPCA  62      -0.211  -1.597  -3.546  0.00  0.00      A    C
ATOM     92  H6X DMPCA  62       0.261  -0.589  -3.502  0.00  0.00      A    H
ATOM     93  H6Y DMPCA  62      -1.308  -1.451  -3.666  0.00  0.00      A    H
ATOM     94  C37 DMPCA  62       0.321  -2.349  -4.774  0.00  0.00      A    C
ATOM     95  H7X DMPCA  62      -0.135  -3.364  -4.802  0.00  0.00      A    H
ATOM     96  H7Y DMPCA  62       1.423  -2.475  -4.681  0.00  0.00      A    H
ATOM     97  C38 DMPCA  62      -0.000  -1.620  -6.084  0.00  0.00      A    C
ATOM     98  H8X DMPCA  62       0.477  -0.613  -6.069  0.00  0.00      A    H
ATOM     99  H8Y DMPCA  62      -1.101  -1.469  -6.152  0.00  0.00      A    H
ATOM    100  C39 DMPCA  62       0.472  -2.383  -7.328  0.00  0.00      A    C
ATOM    101  H9X DMPCA  62       0.008  -3.396  -7.330  0.00  0.00      A    H
ATOM    102  H9Y DMPCA  62       1.576  -2.517  -7.280  0.00  0.00      A    H
ATOM    103 C310 DMPCA  62       0.103  -1.660  -8.630  0.00  0.00      A    C
ATOM    104 H10X DMPCA  62       0.576  -0.651  -8.632  0.00  0.00      A    H
ATOM    105 H10Y DMPCA  62      -1.000  -1.516  -8.663  0.00  0.00      A    H
ATOM    106 C311 DMPCA  62       0.544  -2.427  -9.882  0.00  0.00      A    C
ATOM    107 H11X DMPCA  62       0.078  -3.438  -9.874  0.00  0.00      A    H
ATOM    108 H11Y DMPCA  62       1.649  -2.563  -9.859  0.00  0.00      A    H
ATOM    109 C312 DMPCA  62       0.151  -1.705 -11.178  0.00  0.00      A    C
ATOM    110 H12X DMPCA  62       0.617  -0.693 -11.189  0.00  0.00      A    H
ATOM    111 H12Y DMPCA  62      -0.953  -1.566 -11.198  0.00  0.00      A    H
ATOM    112 C313 DMPCA  62       0.583  -2.470 -12.435  0.00  0.00      A    C
ATOM    113 H13X DMPCA  62       0.119  -3.483 -12.415  0.00  0.00      A    H
ATOM    114 H13Y DMPCA  62       1.688  -2.607 -12.409  0.00  0.00      A    H
ATOM    115 C314 DMPCA  62       0.194  -1.764 -13.738  0.00  0.00      A    C
ATOM    116 H14X DMPCA  62      -0.908  -1.637 -13.795  0.00  0.00      A    H
ATOM    117 H14Y DMPCA  62       0.664  -0.758 -13.790  0.00  0.00      A    H
ATOM    118 H14Z DMPCA  62       0.527  -2.355 -14.616  0.00  0.00      A    H
"
close $lipid_file

# generate the cell.inp file that is included in by the first
# solvated MD simulation namd config file
set fp [open "cell.inp" "w"]
puts $fp "cellbasisvector1 [lindex $basisvec 0] 0 0"
puts $fp "cellbasisvector2 0 [lindex $basisvec 1] 0"
puts $fp "cellbasisvector3 0 0 [lindex $basisvec 2]"
puts $fp "cellorigin $origin"
close $fp
puts "Generated cell.inp."

exit
