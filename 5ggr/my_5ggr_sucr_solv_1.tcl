# VMD script for generating solvated/neutralized system from the protein psf/pdb
#
# cameron f abrams (c) 2018
# cfa22@drexel.edu
# drexel university
# chemical and biological engineering
set csm 0.2
for {set i 0} {$i < $argc} {incr i} {
   if { [lindex $argv $1] == "-cs" } {
     incr i
     set csm [lindex $argv $i]
   }
}
set LOCALFILES {}
set pad 15; # pad in angstroms

set inputname my_5ggr
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

foreach d {0 1 2} {
  lset box 0 $d [expr [lindex $minmax 0 $d] - $pad]
  lset box 1 $d [expr [lindex $minmax 1 $d] + $pad]
  lset basisvec $d [expr [lindex $box 1 $d ] - [lindex $box 0 $d]] 
  lset origin $d [expr 0.5*([lindex $box 1 $d ] + [lindex $box 0 $d])] 
}

$a writepdb prot.pdb
lappend LOCALFILES prot.pdb

set mp [vecsum [$a get mass]]
set z [vecsum [$a get charge]]
set lx [expr [lindex [lindex $box 1] 0] - [lindex [lindex $box 0] 0]]
set ly [expr [lindex [lindex $box 1] 1] - [lindex [lindex $box 0] 1]]
set lz [expr [lindex [lindex $box 1] 2] - [lindex [lindex $box 0] 2]]
set V [expr $lx * $ly * $lz]

set xmin [lindex [lindex $box 0] 0]
set xmax [lindex [lindex $box 1] 0]
set ymin [lindex [lindex $box 0] 1]
set ymax [lindex [lindex $box 1] 1]
set zmin [lindex [lindex $box 0] 2]
set zmax [lindex [lindex $box 1] 2]

puts "$mp $z $lx $ly $lz $V $densamuA3"

# assuming an overall density of 1 g/cc, compute
# the mass fraction (which is assumed equal to the
# volume fraction) of the protein and hence the
# volume of the solvent space
set densgcc 1.0
set pmassfr [expr $mp / (0.6022*$densgcc*$V)]
set Vs [expr (1-$pmassfr)*$V]

# compute number of sucrose molecules
set ns [expr int(6.022e-4 * $csm * $Vs)]
set MWs 163.15

# compute number of water molecules
set MWw 18.0
set nw [expr int( 1.0 / $MWw * (0.6022*$densgcc*$Vs- $MWs*$ns))]

set nna 0
set ncl 0
if { $z > 0 } {
   set ncl [expr round($z)]
} elseif { $z < 0 } {
   set nna [expr round(-($z))]
}

puts "system will have $ns sucrose, $nw waters, $nna sodiums, and $ncl chlorides."

set fp [open "pm-tmp.in" "w"]
puts $fp "
output my_5ggr_sucr.pdb
filetype pdb
tolerance 2.0
structure prot.pdb
  number 1
  center
  fixed 0. 0. 0. 0. 0. 0.
end structure
structure SU.pdb
   resnumbers 2
   number $ns
   inside box $xmin $ymin $zmin $xmax $ymax $zmax
end structure
"
set ws { A B C D E F G I J K M N }; # skip H and L, since those are the protein
set nws [llength $ws]
for {set i 0} {$i < $nws} {incr i} {
  puts $fp "
structure tip3-[lindex $ws $i].pdb
  resnumbers 0
  number [expr int($nw/$nws)]
  inside box $xmin $ymin $zmin $xmax $ymax $zmax
end structure
"
}

if { $nna > 0 } {
   puts $fp "
structure SOD.pdb
   resnumbers 0
   number $nna
   inside box $xmin $ymin $zmin $xmax $ymax $zmax
end structure
"
}

if { $ncl > 0 } {
   puts $fp "
structure CL.pdb
   resnumbers 0
   number $ncl
   inside box $xmin $ymin $zmin $xmax $ymax $zmax
end structure
"
}

for {set i 0} {$i < $nws} {incr i} {
  set water_file [ open tip3-[lindex $ws $i].pdb w ]
  puts $water_file "HETATM    1  H1  TIP3    1       9.626   6.787  12.673  0.00  0.00      W[lindex $ws $i]"
  puts $water_file "HETATM    2  H2  TIP3    1       9.626   8.420  12.673  0.00  0.00      W[lindex $ws $i]"
  puts $water_file "HETATM    3  OH2 TIP3    1      10.203   7.604  12.673  0.00  0.00      W[lindex $ws $i]"
  close $water_file
}

set sod_file [ open SOD.pdb w ]
puts $sod_file  "HETATM    1 SOD  SOD    2        0.000   0.000   0.000  0.00  0.00      I" 
close $sod_file

set cl_file [ open CL.pdb w ] 
puts $cl_file  "HETATM    1 CLA  CLA    2        0.000   0.000   0.000  0.00  0.00      I" 
close $cl_file 

set su_file [ open SU.pdb w ]
puts $su_file \
"ATOM      1  C1  AGLCS   1      -0.409  -2.074  -0.826  0.00  0.00      SU    
ATOM      2  H1  AGLCS   1       0.193  -2.833  -0.283  0.00  0.00      SU    
ATOM      3  O1  AGLCS   1       0.364  -1.285  -1.755  0.00  0.00      SU    
ATOM      4  C5  AGLCS   1      -1.716  -0.165  -0.225  0.00  0.00      SU    
ATOM      5  H5  AGLCS   1      -0.990   0.355  -0.885  0.00  0.00      SU    
ATOM      6  O5  AGLCS   1      -0.926  -1.262   0.240  0.00  0.00      SU    
ATOM      7  C2  AGLCS   1      -1.560  -2.669  -1.550  0.00  0.00      SU    
ATOM      8  H2  AGLCS   1      -2.215  -3.243  -0.860  0.00  0.00      SU    
ATOM      9  O2  AGLCS   1      -1.101  -3.586  -2.551  0.00  0.00      SU    
ATOM     10  HO2 AGLCS   1      -0.538  -3.096  -3.156  0.00  0.00      SU    
ATOM     11  C3  AGLCS   1      -2.325  -1.526  -2.160  0.00  0.00      SU    
ATOM     12  H3  AGLCS   1      -1.621  -0.960  -2.806  0.00  0.00      SU    
ATOM     13  O3  AGLCS   1      -3.408  -1.966  -3.009  0.00  0.00      SU    
ATOM     14  HO3 AGLCS   1      -3.920  -2.651  -2.573  0.00  0.00      SU    
ATOM     15  C4  AGLCS   1      -2.870  -0.650  -0.997  0.00  0.00      SU    
ATOM     16  H4  AGLCS   1      -3.565  -1.295  -0.418  0.00  0.00      SU    
ATOM     17  O4  AGLCS   1      -3.508   0.417  -1.549  0.00  0.00      SU    
ATOM     18  HO4 AGLCS   1      -4.088   0.112  -2.251  0.00  0.00      SU    
ATOM     19  C6  AGLCS   1      -2.106   0.743   0.996  0.00  0.00      SU    
ATOM     20  H61 AGLCS   1      -1.298   1.470   1.223  0.00  0.00      SU    
ATOM     21  H62 AGLCS   1      -2.987   1.353   0.699  0.00  0.00      SU    
ATOM     22  O6  AGLCS   1      -2.512  -0.031   2.077  0.00  0.00      SU    
ATOM     23  HO6 AGLCS   1      -1.759  -0.249   2.632  0.00  0.00      SU    
ATOM     24  O5  BFRUS   2       1.440  -0.967  -2.899  0.00  0.00      SU    
ATOM     25  C2  BFRUS   2       1.765  -1.722  -1.714  0.00  0.00      SU    
ATOM     26  C5  BFRUS   2       1.898  -1.616  -4.067  0.00  0.00      SU    
ATOM     27  H5  BFRUS   2       2.907  -1.212  -4.299  0.00  0.00      SU    
ATOM     28  C6  BFRUS   2       1.011  -1.368  -5.304  0.00  0.00      SU    
ATOM     29  H61 BFRUS   2       1.109  -0.279  -5.499  0.00  0.00      SU    
ATOM     30  H62 BFRUS   2       1.377  -1.936  -6.187  0.00  0.00      SU    
ATOM     31  O6  BFRUS   2      -0.339  -1.633  -4.984  0.00  0.00      SU    
ATOM     32  HO6 BFRUS   2      -0.940  -1.104  -5.514  0.00  0.00      SU    
ATOM     33  C1  BFRUS   2       2.455  -0.693  -0.725  0.00  0.00      SU    
ATOM     34  H11 BFRUS   2       1.704   0.019  -0.321  0.00  0.00      SU    
ATOM     35  H12 BFRUS   2       3.214  -0.062  -1.235  0.00  0.00      SU    
ATOM     36  O1  BFRUS   2       3.126  -1.348   0.377  0.00  0.00      SU    
ATOM     37  HO1 BFRUS   2       2.516  -1.888   0.885  0.00  0.00      SU    
ATOM     38  C3  BFRUS   2       2.772  -2.797  -2.202  0.00  0.00      SU    
ATOM     39  H3  BFRUS   2       3.802  -2.411  -2.354  0.00  0.00      SU    
ATOM     40  O3  BFRUS   2       2.879  -4.040  -1.495  0.00  0.00      SU    
ATOM     41  HO3 BFRUS   2       3.083  -3.772  -0.597  0.00  0.00      SU    
ATOM     42  C4  BFRUS   2       2.149  -3.016  -3.644  0.00  0.00      SU    
ATOM     43  H4  BFRUS   2       1.213  -3.613  -3.689  0.00  0.00      SU    
ATOM     44  O4  BFRUS   2       3.066  -3.716  -4.404  0.00  0.00      SU    
ATOM     45  HO4 BFRUS   2       3.454  -4.420  -3.879  0.00  0.00      SU"    
close $su_file

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
# now, the do_sucr.sh script will run packmol, since it is 
# apparently impossible to invoke from within a TcL script
# after that, the next stage of the psfgen can continue
