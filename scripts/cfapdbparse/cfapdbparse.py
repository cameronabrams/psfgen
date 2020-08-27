import sys
import operator
from datetime import date 
from molecule import Molecule
from cleavage import Cleavage
from mutation import Mutation
from graft import Graft
import argparse
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
            if l.terminated and len(l.residues)>1:
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
    parser=argparse.ArgumentParser()
    print('### cfapdbparser {}'.format(date.today()))
    i=1
    Molecules=[]
    Mut=[]
    Clv=[]
    psfgen='mkpsf.tcl'
    Topo=['top_all36_prot.rtf','top_all36_carb_namd_cfa.rtf','stream/carb/toppar_all36_carb_glycopeptide.str','toppar_water_ions_namd_nonbfixes.str']
    PostMod={}
    PostMod['center_protein']=True
    prefix='x01_'
    fixConflicts=True
    PostMod['do_loop_mc']=False

    parser.add_argument('pdb',nargs='+',metavar='<?.pdb>',type=Molecule,help='name(s) of pdb file to parse; first is treated as the base molecule; others are ')
    parser.add_argument('-topo',metavar='<name>',action='append',default=[],help='additional CHARMM topology files')
    parser.add_argument('-prefix',metavar='<str>',default='x01_',help='output PDB/PSF prefix; each file name will have the format <prefix><pdbcode>.pdb/psf, where <pdbcode> is the 4-letter PDB code of the base molecule.')
    parser.add_argument('-psfgen',metavar='<name>',default='mkpsf.tcl',help='name of TcL script generated as input to VMD/psfgen')
    parser.add_argument('-mut',metavar='X_Y###Z',action='append',default=[],type=Mutation,help='specify mutation.  Format: X is chainID, Y is one-letter residue code to mutate FROM, ### is sequence number (can be any number of digits), and Z is one-letter residue code to mutate TO.  Multiple -mut\'s can be specified.')
    parser.add_argument('-clv',metavar='X###Y',action='append',default=[],type=Cleavage,help='specify cleavage site.  Format: X is parent chain ID, ### is residue number immediately N-terminal to the cleavage site, and Y is the daughter chain ID that will begin immediately C-terminal to cleavage site. Multiple -clv\'s can be specified, each with its own -clv key.')
    parser.add_argument('-gra',metavar='<str>,A:XXX-YYY,ZZZ,C:BBB',action='append',default=[],type=Graft,help='graft resids XXX-YYY of chain A in pdb <str> to chain C of base molecule by overlapping resid ZZZ of chain A of graft and resid BBB of chain C of base')
    parser.add_argument('-graftfile',metavar='<name>',default='',help='input file listing all grafts (as an alternative to issuing multiple -gra arguments)')
    # booleans
    parser.add_argument('-rmi',action='store_true',help='asks psfgen to use the loopMC module to relax modeled-in loops of residues missing from PDB')
    parser.add_argument('-kc',action='store_true',help='ignores SEQADV records indicating conflicts; if unset, residues in conflict are mutated to their proper identities')
    parser.add_argument('-noc',action='store_true',help='do not center the protein at the origin of the coordinate system')
    parser.add_argument('-ror',default='None,None',metavar='<atomselect string>,<atomselect string>',help='two comma-separated, single-quoted atomselect strings to define two groups of atoms whose centers of mass are aligned against the global z-axis')
    parser.add_argument('-v',action='store_true',help='print verbose output during parsing')

    args=parser.parse_args()

    Molecules=args.pdb
    if args.v:
        for m in Molecules:
            m.show()
    Mut=args.mut
    if len(args.topo)>0:
        Topo.extend(args.topo)
    Clv=args.clv
    Gra=args.gra
    if args.graftfile!='':
        with open(args.graftfile,'r') as f:
           for l in f:
               if l[0]!='#':
                   Gra.append(Graft(l))
    prefix=args.prefix
    PostMod['do_loop_mc']=args.rmi
    fixConflicts=~(args.kc)
    psfgen=args.psfgen
    PostMod['center_protein']=~(args.noc)
    if args.ror!='None,None':
        PostMod['reorient_protein']=True
        PostMod['center_protein']=True
        PostMod['reorselstr']=args.ror.split(',')

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
    Loops=Base.writepsfgeninput(psfgen_fp,topologies=Topo,userMutations=Mut,fixConflicts=fixConflicts,prefix=prefix,userGrafts=Gra)
    
    post_pdb=WritePostMods(psfgen_fp,Base.psf_outfile,Base.pdb_outfile,PostMod,Loops)

    Base.Tcl_PrependHeaderToPDB(post_pdb,psfgen_fp)

    psfgen_fp.write('exit\n')
    psfgen_fp.write('### thank you for using cfapdbparser!\n')
    print('### next: vmd -dispdev text -e {}'.format(psfgen))
    print('### will generate {} and {}'.format(Base.psf_outfile,post_pdb))
    
