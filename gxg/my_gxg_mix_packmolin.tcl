# VMD script for generating solvent mixture
#
# cameron f abrams (c) 2019
# cfa22@drexel.edu
# drexel university
# chemical and biological engineering

# 200 mM peptide concentration
set pm 0.2
# 55 wt % ethanol in solvent (other solvent component is water)
set we 0.55
# 30-angstrom cubic box
set L 50.0
set eoh_pdb "my_eoh_q.pdb"
set gxg_pdb "my_gxg_q.pdb"

set seed [exec date +%s]

for {set i 0} {$i < $argc} {incr i} {
   if { [lindex $argv $i] == "-pm" } {
     incr i
     set pm [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-we" } {
     incr i
     set we [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-L" } {
     incr i
     set L [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-seed" } {
     incr i
     set seed [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-eoh_pdb" } {
     incr i
     set eoh_pdb [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-z" } {
     incr i
     set z [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-gxg_pdb" } {
     incr i
     set gxg_pdb [lindex $argv $i]
   }
}
set LOCALFILES {}

if { ! [file exists $eoh_pdb] } {
   puts "ethanol molecule pdb $eoh_pdb not found."
   exit
}

if { ! [file exists $gxg_pdb] } {
   puts "tripeptide molecule pdb $gxg_pdb not found."
   exit
}

set box { { ? ? ? } { ? ? ? } }
set basisvec { ? ? ? }
set origin { ? ? ? }

set hL [expr $L / 2]

foreach d {0 1 2} {
  lset box 0 $d [expr -1*$hL]
  lset box 1 $d $hL
  lset basisvec $d [expr [lindex $box 1 $d ] - [lindex $box 0 $d]] 
  lset origin $d [expr 0.5*([lindex $box 1 $d ] + [lindex $box 0 $d])] 
}

# cubic box for now...
set lx $L
set ly $L
set lz $L

set V [expr $lx * $ly * $lz]

set xmin [lindex [lindex $box 0] 0]
set xmax [lindex [lindex $box 1] 0]
set ymin [lindex [lindex $box 0] 1]
set ymax [lindex [lindex $box 1] 1]
set zmin [lindex [lindex $box 0] 2]
set zmax [lindex [lindex $box 1] 2]

# might have to lower this
set densgcc 0.7

# compute number of GXG molecules
set ngxg [expr int(6.022e-4 * $pm * $V)]

# compute number of counterions
set nna 0
set ncl 0
if { $z > 0 } {
   set ncl [expr round($z)*$ngxg]
} elseif { $z < 0 } {
   set nna [expr round(-$z)*$ngxg]
}

# compute approximate numbers of water and ethanol molecules
set MWw 18.0
set MWe 46.0
set wfw [expr 1.0 - $we]
set molfw [expr $wfw/$MWw/($wfw/$MWw+$we/$MWe)]
set nwifpure [expr int( 1.0 / $MWw * (0.6022*$densgcc*$V))]
set nw [expr int($molfw*$nwifpure)]
set ne [expr int((1-$molfw)*$nwifpure)]


set check_wfw [expr ($nw*$MWw)/($ne*$MWe + $nw*$MWw)]

puts "system will have $ngxg gxg, $nw waters, $ne ethanols, $nna sodium, $ncl chlorides."

set fp [open "pm-tmp.in" "w"]
puts $fp "
output my_gxg_mix.pdb
filetype pdb
seed $seed
tolerance 2.0
structure $gxg_pdb
   resnumbers 1
   number $ngxg
   inside box $xmin $ymin $zmin $xmax $ymax $zmax
end structure
structure $eoh_pdb
   resnumbers 0
   number $ne
   inside box $xmin $ymin $zmin $xmax $ymax $zmax
end structure
"
set ws { A B C D F G I J K M N }
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
close $fp

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
# now, the do_mix.sh script will run packmol, since it is 
# apparently impossible to invoke from within a TcL script
# after that, the next stage of the psfgen can continue
