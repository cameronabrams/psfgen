        fp.write('foreach l $loops {\n')
        fp.write('   set chain [lindex $l 0]\n')
        fp.write('   puts "Relaxing loop $l ($loopindex out of [expr $nloops-1])..."\n')
        fp.write('   set upsel [atomselect $molid "protein and chain $chain and resid [lindex $l 1]"]\n')
        fp.write('   set bnds [$upsel getbonds]\n')
        fp.write('   set an [$upsel get name]\n')
        fp.write('   set ai [$upsel get index]\n')
        fp.write('   set upc -1\n')
        fp.write('   foreach n $an i $ai bl $bnds {\n')
        fp.write('     if { $n == "N" } {\n')
        fp.write('       foreach b $bl {\n')
        fp.write('         if { [lsearch $ai $b] == -1 } {\n')
        fp.write('              set upc $b\n')
        fp.write('         }\n')
        fp.write('       }\n')
        fp.write('     }\n')
        fp.write('   }\n')
        fp.write('   if { $upc == -1 } {\n')
        fp.write('     set firstres [lindex $l 1]\n')
        fp.write('   } else {\n')
        fp.write('     set upres [[atomselect $molid "index $upc"] get resid]\n')
        fp.write('     set firstres $upres\n')
        fp.write('     puts "setting firstres to $upres"\n')
        fp.write('   }\n')
        fp.write('   set msel [atomselect $molid "protein and chain $chain and resid $firstres to [lindex $l 2] and not (resid [lindex $l 2] and name C O)"]\n')
        fp.write('   set atomind [dict create]\n') 
        fp.write('   dict set atomind fa  [[atomselect $molid "protein and chain $chain and resid $firstres and name CA"] get index]\n')
        fp.write('   dict set atomind ca [[atomselect $molid "protein and chain $chain and resid [lindex $l 2] and name CA"] get index]\n')
        fp.write('   dict set atomind c [[atomselect $molid "protein and chain $chain and resid [lindex $l 2] and name C"] get index]\n')
        fp.write('   do_flex_mc $molid $msel $bg atomind mcp [irand_dom 1000 9999] $logid {} {}\n'.format(logevery,logsaveevery))
        fp.write('   set loopindex [expr $loopindex + 1]\n')
        fp.write('}\n')

    if 'do_gly_mc' in PostMod and PostMod['do_gly_mc']:
        nc=1000
        rcut=4.0
        sigma=1.8
        epsilon=0.5
        mctemperature=3.0
        sstop=2.0
        maxanglestep=60.0 # degrees
        if 'gly_mc_params' in PostMod:
            p=PostMod['gly_mc_params']
            nc=nc if 'maxcycles' not in p else p['maxcycles']
            rcut=rcut if 'rcut' not in p else p['rcut']
            sigma=sigma if 'sigma' not in p else p['sigma']
            epsilon=epsilon if 'epsilon' not in p else p['epsilon']
            mctemperature=mctemperature if 'temperature' not in p else p['temperature']
            sstop=sstop if 'sstop' not in p else p['sstop']
            maxanglestep=maxanglestep if 'maxanglestep' not in p else p['maxanglestep']
        fp.write('set mcp [dict create]\n')
        fp.write('dict set mcp nc {}\n'.format(nc))
        fp.write('dict set mcp rcut {}\n'.format(rcut))
        fp.write('dict set mcp sigma {}\n'.format(sigma))
        fp.write('dict set mcp epsilon {}\n'.format(epsilon))
        fp.write('dict set mcp temperature {}\n'.format(mctemperature))
        fp.write('dict set mcp sstop {}\n'.format(sstop))
        fp.write('dict set mcp dstop -1\n')
        fp.write('dict set mcp maxanglestep {}\n'.format(maxanglestep))
        fp.write('set bg [atomselect $molid "noh"]\n')
        fp.write('set glycan_segs [list '+' '.join(GlycanSegs)+']\n')
        fp.write('set ng [llength $glycan_segs]\n')
        fp.write('set gi 1\n')
        fp.write('foreach g $glycan_segs {\n')
        fp.write('   set sel [atomselect $molid "segname $g"]\n')
        fp.write('   set rid [$sel get resid]\n')
        fp.write('   set root [lindex [lsort -unique -real $rid] 0]\n')
        fp.write('   set atomind [dict create]\n') 
        fp.write('   dict set atomind fa  [[atomselect $molid "segname $g and name C1 and resid $root"] get index]\n')
        fp.write('   dict set atomind ca  -1\n')
        fp.write('   dict set atomind c  -1\n')
        fp.write('   puts "Relaxing glycan $g ($gi/$ng) rootres $root..."\n')
        fp.write('   do_flex_mc $molid $msel $bg atomind mcp [irand_dom 1000 9999] $logid {} {}\n'.format(logevery,logsaveevery))
        fp.write('   set gi [expr $gi + 1]\n')
        fp.write('}\n')