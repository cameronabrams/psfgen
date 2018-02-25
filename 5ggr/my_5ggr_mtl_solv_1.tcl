# VMD script for generating solvated/neutralized system from the protein psf/pdb
#
# cameron f abrams (c) 2018
# cfa22@drexel.edu
# drexel university
# chemical and biological engineering
set cmm 0.2
for {set i 0} {$i < $argc} {incr i} {
   if { [lindex $argv $1] == "-cm" } {
     incr i
     set cmm [lindex $argv $i]
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

set seed [exec date +%s]

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

puts "$mp $z $lx $ly $lz $V"

# assuming an overall density of 1 g/cc, compute
# the mass fraction (which is assumed equal to the
# volume fraction) of the protein and hence the
# volume of the solvent space
set densgcc 1.0
set pmassfr [expr $mp / (0.6022*$densgcc*$V)]
set Vs [expr (1-$pmassfr)*$V]

# compute number of mannitol molecules
set nm [expr int(6.022e-4 * $cmm * $Vs)]
set MWm 182.172

# compute number of water molecules
set MWw 18.0
set nw [expr int( 1.0 / $MWw * (0.6022*$densgcc*$Vs- $MWm*$nm))]

set nna 0
set ncl 0
if { $z > 0 } {
   set ncl [expr round($z)]
} elseif { $z < 0 } {
   set nna [expr round(-($z))]
}

puts "system will have $nm mannitol, $nw waters, $nna sodiums, and $ncl chlorides."

set fp [open "pm-tmp.in" "w"]
puts $fp "
output my_5ggr_mtl.pdb
filetype pdb
seed $seed
tolerance 2.0
structure prot.pdb
  number 1
  center
  fixed 0. 0. 0. 0. 0. 0.
end structure
structure MTL.pdb
   resnumbers 0
   number $nm
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

set mtl_file [ open MTL.pdb w ]
puts $mtl_file \
"ATOM      1  C1  DMANM   1      -0.023  -2.035   2.486  0.00  0.00      M
ATOM      2  H11 DMANM   1      -0.007  -3.055   2.060  0.00  0.00      M
ATOM      3  H12 DMANM   1       1.091  -1.855   2.795  0.00  0.00      M
ATOM      4  O1  DMANM   1      -0.891  -2.150   3.642  0.00  0.00      M
ATOM      5  HO1 DMANM   1      -1.500  -2.828   3.490  0.00  0.00      M
ATOM      6  C2  DMANM   1      -0.556  -0.952   1.599  0.00  0.00      M
ATOM      7  H2  DMANM   1      -1.571  -1.206   1.397  0.00  0.00      M
ATOM      8  O2  DMANM   1      -0.662   0.328   2.300  0.00  0.00      M
ATOM      9  HO2 DMANM   1      -1.642   0.375   2.439  0.00  0.00      M
ATOM     10  C3  DMANM   1       0.263  -0.675   0.327  0.00  0.00      M
ATOM     11  H3  DMANM   1       1.345  -0.626   0.600  0.00  0.00      M
ATOM     12  O3  DMANM   1       0.085  -1.808  -0.568  0.00  0.00      M
ATOM     13  HO3 DMANM   1      -0.352  -2.496  -0.035  0.00  0.00      M
ATOM     14  C4  DMANM   1      -0.092   0.642  -0.409  0.00  0.00      M
ATOM     15  H4  DMANM   1       0.126   1.468   0.316  0.00  0.00      M
ATOM     16  O4  DMANM   1      -1.495   0.698  -0.743  0.00  0.00      M
ATOM     17  HO4 DMANM   1      -1.889   1.155   0.021  0.00  0.00      M
ATOM     18  C5  DMANM   1       0.709   0.873  -1.727  0.00  0.00      M
ATOM     19  H5  DMANM   1       0.424   0.081  -2.506  0.00  0.00      M
ATOM     20  O5  DMANM   1       2.111   0.841  -1.398  0.00  0.00      M
ATOM     21  HO5 DMANM   1       2.500   0.114  -1.882  0.00  0.00      M
ATOM     22  C6  DMANM   1       0.493   2.272  -2.320  0.00  0.00      M
ATOM     23  H61 DMANM   1      -0.594   2.513  -2.275  0.00  0.00      M
ATOM     24  H62 DMANM   1       1.172   3.015  -1.826  0.00  0.00      M
ATOM     25  O6  DMANM   1       0.753   2.277  -3.730  0.00  0.00      M
ATOM     26  HO6 DMANM   1       0.203   3.034  -4.053  0.00  0.00      M"    
close $mtl_file

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
