mol new my_5u7o-davei-l7_i.psf
#mol addfile my_5u7o-davei-l7_i.pdb
mol addfile sol-stage3.dcd waitfor all
mol addfile sol-stage3_r1.dcd waitfor all

set REMAKE_DX 0

set ref [atomselect top "protein and chain G K R A B C"]

$ref frame 0
set a [atomselect top "protein and chain G K R A B C"]
for { set i 1 } { $i < [molinfo top get numframes] } { incr i } {
  $a frame $i
  $a move [measure fit $a $ref]
}

foreach c {G K R A B C} {
  if { ! [file exists map_${c}.dx] || $REMAKE_DX } {
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
mol selection {chain K}
mol material AOChalky
mol addrep top
mol representation Isosurface $isoden 2 0 0 1 1
mol color ColorID 6
mol selection {chain R}
mol material AOChalky
mol addrep top
mol representation Isosurface $isoden 3 0 0 1 1
mol color ColorID 8
mol selection {chain A}
mol material AOChalky
mol addrep top
mol representation Isosurface $isoden 4 0 0 1 1
mol color ColorID 8
mol selection {chain B}
mol material AOChalky
mol addrep top
mol representation Isosurface $isoden 5 0 0 1 1
mol color ColorID 8
mol selection {chain C}
mol material Transparent
mol addrep top
mol representation Tube 0.3 12.0
mol color ColorID 6
mol selection {chain G}
mol material Opaque
mol addrep top
mol representation Tube 0.3 12.0
mol color ColorID 8
mol selection {chain C}
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

set lf 250
mol drawframes top 6 [list $lf]
mol drawframes top 7 [list $lf]
mol drawframes top 8 [list $lf]
mol drawframes top 9 [list $lf]
mol drawframes top 10 [list $lf]
mol drawframes top 11 [list $lf]

set vp {{{1 0 0 20.4943} {0 1 0 -5.76293} {0 0 1 -0.409865} {0 0 0 1}} {{1.0 0.0 0.0 0} {0.0 0.0 -1.0 0} {0.0 1.0 0.0 0} {0 0 0 1}} {{0.021842 0 0 0} {0 0.021842 0 0} {0 0 0.021842 0} {0 0 0 1}} {{1 0 0 -0.47} {0 1 0 0.19} {0 0 1 0.44} {0 0 0 1}}}
molinfo top set {center_matrix rotate_matrix scale_matrix global_matrix} $vp
display projection orthographic
display resize 1000 1000
color Display {Background} white
axes location off

render TachyonInternal scene.tga

exit 
