mol new my_5vn3_i.psf
#mol addfile my_5u7o-davei-l7_i.pdb
mol addfile sol-stage2.dcd waitfor all

set REMAKE_DX 0

set ref [atomselect top "protein and chain G I J B A D"]

$ref frame 0
set a [atomselect top "protein and chain G I J B A D"]
for { set i 1 } { $i < [molinfo top get numframes] } { incr i } {
  $a frame $i
  $a move [measure fit $a $ref]
}

foreach c {G I J B A D} {
  if { $REMAKE_DX || ! [file exists map_${c}.dx] } {
    volmap density [atomselect top "chain $c ${c}S"] -allframes -res 2.0 -radscale 9.0 -mol top -o map_${c}.dx
  } else {
    mol addfile map_${c}.dx
  }
}

mol delrep 0 top

set isoden 0.04

mol representation Isosurface $isoden 0 0 0 1 1
mol color ColorID 6
mol selection {chain G}
mol material Transparent
mol addrep top
mol representation Isosurface $isoden 1 0 0 1 1
mol color ColorID 6
mol selection {chain I}
mol material AOChalky
mol addrep top
mol representation Isosurface $isoden 2 0 0 1 1
mol color ColorID 6
mol selection {chain J}
mol material AOChalky
mol addrep top
mol representation Isosurface $isoden 3 0 0 1 1
mol color ColorID 8
mol selection {chain B}
mol material AOChalky
mol addrep top
mol representation Isosurface $isoden 4 0 0 1 1
mol color ColorID 8
mol selection {chain A}
mol material AOChalky
mol addrep top
mol representation Isosurface $isoden 5 0 0 1 1
mol color ColorID 8
mol selection {chain D}
mol material Transparent
mol addrep top
mol representation Tube 0.3 12.0
mol color ColorID 6
mol selection {chain G}
mol material Opaque
mol addrep top
mol representation Tube 0.3 12.0
mol color ColorID 8
mol selection {chain D}
mol material Opaque
mol addrep top

mol representation Licorice 0.5 15.0 12.0
mol color ColorID 7
mol selection {segname X or (segname X1 and within 2.0 of segname X)}
mol material Goodsell
mol addrep top
mol representation Licorice 0.3 15.0 12.0
mol color ColorID 16
mol selection {segname X1 or (segname X2 and within 2.0 of segname X1)}
mol material Goodsell
mol addrep top
mol representation Licorice 0.3 15.0 12.0
mol color ColorID 0
mol selection {segname X2 or (segname XT and within 2.0 of segname X2)}
mol material Goodsell
mol addrep top
mol representation Licorice 0.3 15.0 12.0
mol color ColorID 1
mol selection {segname XT}
mol material Goodsell
mol addrep top

set lf 0
mol drawframes top 6 [list $lf]
mol drawframes top 7 [list $lf]
mol drawframes top 8 [list $lf]
mol drawframes top 9 [list $lf]
mol drawframes top 10 [list $lf]
mol drawframes top 11 [list $lf]
set vp {{{1 0 0 -7.06653} {0 1 0 -3.44793} {0 0 1 4.19107} {0 0 0 1}} {{-0.788031 0.615315 0.0198297 0} {-0.00644338 -0.0404521 0.999161 0} {0.615599 0.787241 0.0358427 0} {0 0 0 1}} {{0.021842 0 0 0} {0 0.021842 0 0} {0 0 0.021842 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}
molinfo top set {center_matrix rotate_matrix scale_matrix global_matrix} $vp
display projection orthographic
display resize 1000 1000
color Display {Background} white
axes location off

render TachyonInternal scene.tga

exit
