import sys
import operator
from datetime import date 
from molecule import Molecule
from cleavage import Cleavage
from mutation import Mutation
''' 
    Parses experimental PDB to build input file for VMD/psfgen
    Cameron F Abrams
    cfa22@drexel.edu
'''
class PostMod:
    def __init__(self,center_protein=False,reorient_protein=False,reorselstr='',do_loop_mc=False,Loops=[]):
        self.ctr=center_protein
        return self

def WritePostMods(fp,psf,pdb,PostMod,Loops):
    prefix=pdb[:]
    prefix=prefix.replace('.pdb','')
    fp.write('### Post modifications follow:\n')
    fp.write('mol delete top\n')
    fp.write('mol new {}\n'.format(psf))
    fp.write('set molid [molinfo top get id]\n')
    fp.write('mol addfile {}\n'.format(pdb))
    ctr=False
    if 'center_protein' in PostMod:
       ctr=PostMod['center_protein']
    if ctr:
        fp.write('set a [atomselect top "all"]\n')
        fp.write('set or [measure center $a weight mass]\n')
        fp.write('$a moveby [vecscale -1 $or]\n')
        reor=False
        if 'reorient_protein' in PostMod:
            reor=PostMod['reorient_protein']
        if reor:
            fp.write('set ca [measure center [atomselect top "protein and {}"] weight mass]\n'.format(PostMod['reorselstr'][0]))
            fp.write('set cb [measure center [atomselect top "protein and {}"] weight mass]\n'.format(PostMod['reorselstr'][1]))
            fp.write('set pi 3.415928\n')
            fp.write('set dv [vecsub $ca $cb]\n')
            fp.write('set d [veclength $dv]\n')
            fp.write('set cp [expr [lindex $dv 0]/$d]\n')
            fp.write('set sp [expr [lindex $dv 1]/$d]\n')
            fp.write('set p [expr acos($cp)]\n')
            fp.write('if {[expr $sp < 0.0]} {\n')
            fp.write('  set p [expr 2*$pi-$p]\n')
            fp.write('}\n')
            fp.write('set ct [expr [lindex $dv 2]/$d]\n')
            fp.write('set t [expr acos($ct)]\n')
            fp.write('$a move [transaxis z [expr -1 * $p] rad]\n')
            fp.write('$a move [transaxis y [expr -1 * $t] rad]\n')
    dlmc=False
    if 'do_loop_mc' in PostMod:
        if PostMod['do_loop_mc']:
           dlmc=True
    if dlmc:
        fp.write('set loops {\n')
        for l in Loops:
            if l.terminated:
                fp.write('{{ {} {} {} }}\n'.format(l.chainID,l.residues[0].resseqnum,l.residues[-1].resseqnum))
        fp.write('           }\n')
        # create loops list { { }, { }, ...}
        fp.write('set nc 1000\n')
        fp.write('set rcut 3.0\n')
        fp.write('set r0 1.5\n')
        fp.write('set temperature 3.0\n')
        fp.write('set k 10.0\n')
        fp.write('set bg [atomselect $molid "noh"]\n')
        fp.write('set loopindex 0\n')
        fp.write('set nloops [llength $loops]\n')
        fp.write('foreach l $loops {\n')
        fp.write('   set chain [lindex $l 0]\n')
        fp.write('   puts "Relaxing loop $loopindex out of $nloops"\n')
        fp.write('   set residueList [[atomselect $molid "chain $chain and resid [lindex $l 1] to [lindex $l 2] and name CA"] get residue]\n')
        fp.write('   do_loop_mc $residueList $chain $molid $k $r0 $bg $rcut $nc $temperature [irand_dom 1000 9999] $logid\n')
        fp.write('   set loopindex [expr $loopindex + 1]\n')
        fp.write('}\n')

    new_pdb_out=prefix+'_mod.pdb'
    fp.write('$a writepdb {}\n'.format(new_pdb_out))
    return new_pdb_out

if __name__=='__main__':

    print('### cfapdbparser {}'.format(date.today()))
    i=1
    Molecules=[]
    Mut=[]
    Clv=[]
    psfgen='mkpsf.tcl'
    topo=['top_all36_prot.rtf','top_all36_carb_namd_cfa.rtf','stream/carb/toppar_all36_carb_glycopeptide.str','toppar_water_ions_namd_nonbfixes.str']
    PostMod={}
    PostMod['center_protein']=True
    while i<len(sys.argv):
        if sys.argv[i]=='-pdb':
            i+=1
            j=0
            while i<len(sys.argv) and sys.argv[i][0]!='-':
                m=Molecule(j,sys.argv[i])
                if m!=0:
                   print('###',m)
                   Molecules.append(m)
                   j+=1
                i+=1
            if i<len(sys.argv) and sys.argv[i][0]=='-':
                i-=1
        elif sys.argv[i]=='-mut':
            i+=1
            while i<len(sys.argv) and sys.argv[i][0]!='-':
                Mut.append(Mutation(sys.argv[i]))
                i+=1
            if i<len(sys.argv) and sys.argv[i][0]=='-':
                i-=1
        elif sys.argv[i]=='-cleave':
            i+=1
            while i<len(sys.argv) and sys.argv[i][0]!='-':
                Clv.append(Cleavage(sys.argv[i]))
                i+=1
            if i<len(sys.argv) and sys.argv[i][0]=='-':
                i-=1
        elif sys.argv[i]=='-top':
            i+=1
            while i<len(sys.argv) and sys.argv[i][0]!='-':
                topo.append(sys.argv[i])
                i+=1
            if i<len(sys.argv) and sys.argv[i][0]=='-':
                i-=1
        elif sys.argv[i]=='-do_loop_mc':
            PostMod['do_loop_mc']=True
        elif sys.argv[i]=='-no_center':
            PostMod['center_protein']=False
        elif sys.argv[i]=='-reorient_protein':
            PostMod['reorient_protein']=True
            PostMod['reorselstr']=[]
            i+=1
            while i<len(sys.argv) and sys.argv[i][0]!='-':
                PostMod['reorselstr'].append(sys.argv[i])
                i+=1
            if sys.argv[i][0]=='-':
                i-=1
        elif sys.argv[i]=='-psfgen':
            i+=1
            psfgen=sys.argv[i]
        i+=1

    print('### this will generate the following psfgen script: {}'.format(psfgen))

    psfgen_fp=open(psfgen,'w')
    psfgen_fp.write('### This is an automatically generated psfgen input file\n')
    psfgen_fp.write('### created using cfapdbparser.py on {}\n'.format(date.today()))
    psfgen_fp.write('### as part of the psfgen repository\n')
    psfgen_fp.write('### github.com/cameronabrams/psgen\n')
    psfgen_fp.write('### questions to cfa22@drexel.edu\n')
    psfgen_fp.write('### command: python3 ')
    for a in sys.argv:
        psfgen_fp.write('{} '.format(a))
    psfgen_fp.write('\n')
    
    Base=Molecules[0]
    if len(Clv)>0:
        Base.CleaveChains(Clv)
    Loops=Base.writepsfgeninput(psfgen_fp,Mut,topo)
    
    post_pdb=WritePostMods(psfgen_fp,Base.psf_outfile,Base.pdb_outfile,PostMod,Loops)

    Base.Tcl_PrependHeaderToPDB(post_pdb,psfgen_fp)

    psfgen_fp.write('exit\n')
    psfgen_fp.write('### thank you for using cfapdbparser!\n')
    print('### next: vmd -dispdev text -e {}'.format(psfgen))
    print('### will generate {} and {}'.format(Base.psf_outfile,post_pdb))
    
