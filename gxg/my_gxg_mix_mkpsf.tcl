set firstpatch GLYP
set lastpatch CTER
for {set i 0} {$i < $argc} {incr i} {
   if { [lindex $argv $i] == "-firstpatch" } {
      incr i
      set firstpatch [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-lastpatch" } {
      incr i
      set lastpatch [lindex $argv $i]
   }
}

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_na.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb.rtf
topology $env(HOME)/charmm/toppar/top_all36_lipid.rtf
topology $env(HOME)/charmm/toppar/top_all36_cgenff.rtf
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd.str

mol new my_gxg_mix.pdb ; # output of packmol
set ngxgs [expr [[atomselect top "chain A and name CA"] num]/3]
set ncl [expr [[atomselect top "name CLA"] num]]
set nna [expr [[atomselect top "name SOD"] num]]

puts "$ngxgs gxgs found in my_gxg_mix.pdb"
for {set g 0} {$g < $ngxgs} {incr g} {
   set r0 [expr $g * 3]
   set r1 [expr ($g+1)*3-1]
   set tgxg [atomselect top "chain A and residue $r0 to $r1"]
   $tgxg set segname G${g}
   $tgxg writepdb ${g}.pdb
   segment G${g} {
      first $firstpatch
      last $lastpatch
      pdb ${g}.pdb
   }
   coordpdb ${g}.pdb G${g}
}

if { $ncl > 0 } {
  [atomselect top "name CLA"] writepdb cla.pdb
  segment I2 {
    pdb cla.pdb
  }
  coordpdb cla.pdb I2
}

if { $nna > 0 } {
  [atomselect top "name SOD"] writepdb sod.pdb
  segment I1 {
    pdb sod.pdb
  }
  coordpdb sod.pdb I1
}

# cosolvent
set segs { E WA WB WC WD WF WG WI WJ WK WM WN }
set chns { E W  W  W  W  W  W  W  W  W  W  W }
foreach sn $segs cn $chns {
  set s [atomselect top "segname $sn"] 
  $s set chain $cn
  $s writepdb "${sn}.pdb"
  segment $sn {
    first none
    last none
#    if { [string index $sn 0] == "W" } {
      auto none
#    }
    pdb ${sn}.pdb
  }
  coordpdb ${sn}.pdb $sn
}

regenerate angles dihedrals

set outputname my_gxg_mix_solv
writepsf "${outputname}.psf"
writepdb "${outputname}.pdb"

quit
