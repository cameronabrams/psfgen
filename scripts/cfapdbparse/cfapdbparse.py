"""
    Parses PDB/mmCIF file(s) to build input file for VMD/psfgen
    Cameron F Abrams
    cfa22@drexel.edu

"""

import sys
import operator
import argparse
from datetime import date 
from molecule import Molecule
from cleavage import Cleavage
from mutation import Mutation
from ssbond import SSBond
from graft import Graft
from crot import Crot
from attach import Attach
from link import Link
from atom import _PDBAtomNameDict_
from residue import Residue, _PDBResName123_, _pdb_glycans_, _pdb_ions_, _ResNameDict_PDB_to_CHARMM_, _ResNameDict_CHARMM_to_PDB_, get_residue

def WritePostMods(fp,psf,pdb,PostMod,Loops,GlycanSegs):
    """ Writes TcL/VMD commands that encode modifications once the
        the base psfgen structure has been written.  'PostMods' include
        things like commands to center the protein, relax model-built loops, etc.
    """
    logfile=''
    logevery=1
    if 'log_dcd_file' in PostMod:
        logfile=PostMod['log_dcd_file']
    if 'log_every' in PostMod:
        logevery=PostMod['log_every']
    logdcd=len(logfile)>0

    prefix=pdb[:]
    prefix=prefix.replace('.pdb','')
    fp.write('### Post modifications follow:\n')
    fp.write('mol delete top\n')
    fp.write('mol new {}\n'.format(psf))
    fp.write('set molid [molinfo top get id]\n')
    fp.write('mol addfile {}\n'.format(pdb))
    if logdcd:
        fp.write('### logging enabled\n')
        fp.write('mol new {}\n'.format(psf))
        fp.write('mol addfile {}\n'.format(pdb))
        fp.write('set logid [molinfo top get id]\n')
        fp.write('mol top $molid\n')
    else:
        fp.write('set logid -1\n')
    if 'center_protein' in PostMod and PostMod['center_protein']:
        fp.write('set a [atomselect $molid "all"]\n')
        fp.write('set or [measure center $a weight mass]\n')
        fp.write('$a moveby [vecscale -1 $or]\n')
        if logdcd:
            fp.write('set la [atomselect $logid "all"]\n')
            fp.write('$la moveby [vecscale -1 $or]\n')
        if 'reorient_protein' in PostMod and PostMod['reorient_protein']:
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
            if logdcd:
                fp.write('$la move [transaxis z [expr -1 * $p] rad\n')
                fp.write('$la move [transaxis y [expr -1 & $t] rad\n')
    for crot in PostMod['Crot']:
        fp.write(crot.psfgen_str())
        if logdcd:
            fp.write('log_addframe $molid $logid\n')
    if 'do_loop_mc' in PostMod and PostMod['do_loop_mc']:
        fp.write('set loops {\n')
        for l in Loops:
            if l.term and len(l.residues)>1:
                fp.write('{{ {} {} {} }}\n'.format(l.replica_chainID,l.residues[0].resseqnum,l.residues[-1].resseqnum))
        fp.write('           }\n')
        fp.write('set nc 1000\n')
        fp.write('set rcut 4.0\n')
        fp.write('set sigma 1.8\n')
        fp.write('set epsilon 0.5\n')
        fp.write('set r0 1.5\n')
        fp.write('set temperature 3.0\n')
        fp.write('set k 100.0\n')
        fp.write('set bg [atomselect $molid "noh"]\n')
        fp.write('set loopindex 0\n')
        fp.write('set nloops [llength $loops]\n')
        fp.write('foreach l $loops {\n')
        fp.write('   set chain [lindex $l 0]\n')
        fp.write('   puts "Relaxing loop $loopindex out of $nloops"\n')
        fp.write('   set residueList [[atomselect $molid "chain $chain and resid [lindex $l 1] to [lindex $l 2] and name CA"] get residue]\n')
        fp.write('   do_loop_mc $residueList $chain $molid $k $r0 $bg $sigma $epsilon $rcut $nc $temperature [irand_dom 1000 9999] $logid {}\n'.format(logevery))
        fp.write('   set loopindex [expr $loopindex + 1]\n')
        fp.write('}\n')

    if 'do_gly_mc' in PostMod and PostMod['do_gly_mc']:
        fp.write('set nc 1000\n')
        fp.write('set rcut 4.0\n')
        fp.write('set sigma 1.8\n')
        fp.write('set epsilon 0.5\n')
        fp.write('set temperature 3.0\n')
        fp.write('set bg [atomselect $molid "noh"]\n')
        fp.write('set glycan_segs [list '+' '.join(GlycanSegs)+']\n')
        fp.write('foreach g $glycan_segs {\n')
        fp.write('   set sel [atomselect $molid "segname $g"]\n')
        fp.write('   set rid [$sel get resid]\n')
        fp.write('   set root [lindex [lsort -unique -real $rid] 0]\n')
        fp.write('   set fa [[atomselect $molid "segname $g and name C1 and resid $root"] get index]\n')
        fp.write('   puts "Relaxing glycan $g rootres $root rootatom $fa..."\n')
        fp.write('   do_flex_mc $molid $sel $fa 0 -1 -1 $bg $sigma $epsilon $rcut $nc $temperature [irand_dom 1000 9999] $logid {}\n'.format(logevery))
        fp.write('}\n')
    new_pdb_out=prefix+'_mod.pdb'
    fp.write('$a writepdb {}\n'.format(new_pdb_out))
    if logdcd:
        fp.write('set loga [atomselect $logid all]\n')
        fp.write('animate write dcd {} waitfor all sel $loga $logid\n'.format(logfile))
        fp.write('mol delete $logid\n')
    return new_pdb_out

def WriteHeaders(fp,charmm_topologies,local_topologies):
    fp.write('if {![info exists PSFGEN_BASEDIR]} {\n'+\
	  '    if {[info exists env(PSFGEN_BASEDIR)]} {\n'+\
	  '        set PSFGEN_BASEDIR $env(PSFGEN_BASEDIR)\n'+\
	  '    } else {\n'+\
	  '        set PSFGEN_BASEDIR $env(HOME)/research/psfgen\n'+\
	  '    }\n'+\
	  '}\n'+
          'set LOCAL_TOPPARDIR $PSFGEN_BASEDIR/charmm\n')
    fp.write('if {![info exists CHARMM_TOPPARDIR]} {\n'+\
	  '    if {[info exists env(CHARMM_TOPPARDIR)]} {\n'+\
	  '        set TOPPARDIR $env(CHARMM_TOPPARDIR)\n'+\
	  '    } else {\n'+\
	  '        set TOPPARDIR $env(HOME)/charmm/toppar\n'+\
	  '    }\n'+\
	  '}\n')
    fp.write('source ${PSFGEN_BASEDIR}/src/loopmc.tcl\n')
    fp.write('source ${PSFGEN_BASEDIR}/scripts/vmdrc.tcl\n')
    fp.write('package require psfgen\n')
    for t in charmm_topologies:
        fp.write('topology $TOPPARDIR/{}\n'.format(t))
    for t in local_topologies:
        fp.write('topology $LOCAL_TOPPARDIR/{}\n'.format(t))
    fp.write('pdbalias residue HIS HSD\n')
    fp.write('pdbalias atom ILE CD1 CD\n')
    fp.write('pdbalias residue NAG BGNA\n')
    fp.write('pdbalias atom BGNA C7 C\n')
    fp.write('pdbalias atom BGNA O7 O\n')
    fp.write('pdbalias atom BGNA C8 CT\n')
    fp.write('pdbalias atom BGNA N2 N\n')
    fp.write('pdbalias residue SIA ANE5\n')
    fp.write('pdbalias atom ANE5 C10 C\n')
    fp.write('pdbalias atom ANE5 C11 CT\n')
    fp.write('pdbalias atom ANE5 N5 N\n')
    fp.write('pdbalias atom ANE5 O1A O11\n')
    fp.write('pdbalias atom ANE5 O1B O12\n')
    fp.write('pdbalias atom ANE5 O10 O\n')

    for k,v in _ResNameDict_PDB_to_CHARMM_.items():
        fp.write('set RESDICT({}) {}\n'.format(k,v))
    for k,v in _PDBAtomNameDict_.items():
        fp.write('set ANAMEDICT({}) {}\n'.format(k,v))

def MrgCmdLineAndFileContents(cl_list,filename,typ):
    if filename!='':
        with open(filename,'r') as f:
           for l in f:
               if l[0]!='#':
                   cl_list.append(typ(l))
    return cl_list

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    print('cfapdbparser {} / python {}'.format(date.today(),sys.version.replace('\n',' ').split(' ')[0]))
    i=1
    Molecules=[]
    Mut=[]
    Clv=[]
    Uss=[]
    UIC=[]
    psfgen='mkpsf.tcl'
    CTopo=['top_all36_prot.rtf','stream/carb/toppar_all36_carb_glycopeptide.str']
    LocTopo=['top_all36_carb.rtf','toppar_water_ions.str']
    PostMod={}
    PostMod['center_protein']=True
    prefix='x01_'
    fixConflicts=True
    PostMod['do_loop_mc']=False
    PostMod['do_gly_mc']=False
    PostMod['Crot']=[]

    parser.add_argument('pdbcif',nargs='+',metavar='<?.pdb|cif>',type=str,help='name(s) of pdb or CIF file to parse; first is treated as the base molecule; others are not considered (for now)')
    parser.add_argument('-topo',metavar='<name>',action='append',default=[],help='additional CHARMM topology files')
    parser.add_argument('-prefix',metavar='<str>',default='x01_',help='output PDB/PSF prefix; each file name will have the format <prefix><pdbcode>.pdb/psf, where <pdbcode> is the 4-letter PDB code of the base molecule.')
    parser.add_argument('-psfgen',metavar='<name>',default='mkpsf.tcl',help='name of TcL script generated as input to VMD/psfgen')
    parser.add_argument('-ignore',metavar='X',action='append',default=[],type=str,help='Specify a chain to ignore.  Multiple -ignore switches can be used to ignore more than one chain.')
    parser.add_argument('-mut',metavar='X_Y###Z',action='append',default=[],type=Mutation,help='specify mutation.  Format: X is chainID, Y is one-letter residue code to mutate FROM, ### is sequence number (can be any number of digits), and Z is one-letter residue code to mutate TO.  Multiple -mut\'s can be specified.  Mutations are automatically replicated if there are BIOMT transformations.')
    parser.add_argument('-mutfile',metavar='<name>',default='',help='input file listing all mutations (as an alternative to issuing multiple -mut arguments)')
    parser.add_argument('-clv',metavar='X###Y',action='append',default=[],type=Cleavage,help='specify cleavage site.  Format: X is parent chain ID, ### is residue number immediately N-terminal to the cleavage site, and Y is the daughter chain ID that will begin immediately C-terminal to cleavage site. Multiple -clv\'s can be specified, each with its own -clv key.')
    parser.add_argument('-clvfile',metavar='<name>',default='',help='input file listing all cleavages (as an alternative to issuing multiple -clv arguments)')
    parser.add_argument('-gra',metavar='<str>,A:XXX-YYY,ZZZ,C:BBB',action='append',default=[],type=Graft,help='graft resids XXX-YYY of chain A in pdb <str> to chain C of base molecule by overlapping resid ZZZ of chain A of graft and resid BBB of chain C of base.  Grafts are automatically replicated if there are BIOMT transformations.')
    parser.add_argument('-grafile',metavar='<name>',default='',help='input file listing all grafts (as an alternative to issuing multiple -gra arguments)')
    parser.add_argument('-att',metavar='<str>,A:XXX-YYY,ZZZ,B:QQQ,C:BBB',action='append',default=[],type=Attach,help='attach resids XXX-YYY of chain A using resid ZZZ (between XXX and YYY) in pdb <str> to chain C of base molecule at resid BBB by aligning resid QQQ of chain B from source to resid BBB of chain C of base.  CURRENTLY UNIMPLEMENTED!!!')
    parser.add_argument('-attfile',metavar='<name>',default='',help='input file listing all attachments (as an alternative to issuing multiple -att arguments)')
    parser.add_argument('-crot',metavar='<str>,A,XXX[,YYY],###',default=[],action='append',type=Crot,help='specify rotation about a specific torsion angle.  <str> is one of phi, psi, omega, chi1, or chi2.  A is the chainID, XXX is the resid of owner of torson, and YYY (if given) marks the end of the sequence C-terminal to XXX that is reoriented by a backbone rotation. ### is the degrees of rotation.  C-rotations are automatically replicated if there are BIOMT transformations.')
    parser.add_argument('-crotfile',metavar='<name>',default='',help='input file listing all torsion rotations requested (as an alternative to issuing multiple -crot arguments)')
    parser.add_argument('-ssbond',metavar='X_###-Y_###',default=[],action='append',type=SSBond,help='Specify a disulfide bond not in the PDB file; X,Y are chainIDs and ### are resids; if residues are not CYS in wt or by mutations, there is no effect.  Because SSBonds can join chains together, they are NOT automatially replicated if there are BIOMT transformations.')
    parser.add_argument('-ssfile',metavar='<name>',default='',help='input file listing all disulfide bonds to add that are not already in the PDB file (as an alternative to issuing multiple -ssbond arguments)')
    parser.add_argument('-link',metavar='string',default=[],action='append',type=Link,help='PDB-format LINK record; must have exact spacing')
    parser.add_argument('-linkfile',metavar='<name>',default='',help='input file with PDB-format LINK records the user would like to enforce that are not in the RCSB PDB file')
    parser.add_argument('-logdcd',metavar='<name>.dcd',default='',help='name of dcd logging file')
    parser.add_argument('-logevery',metavar='<int>',default=1,help='number of MC accepts between successive frame logging')
    # booleans
    parser.add_argument('-rmi',action='store_true',help='asks psfgen to use the loopMC module to relax modeled-in loops of residues missing from PDB')
    parser.add_argument('-grel',action='store_true',help='asks psfgen to use the loopMC module to relax modeled-in glycans missing from PDB')
    parser.add_argument('-kc',action='store_true',help='ignores SEQADV records indicating conflicts; if unset, residues in conflict are mutated to their proper identities')
    parser.add_argument('-rem',action='store_true',help='revert engineered mutations listed in SEQADV records')
    parser.add_argument('-noc',action='store_true',help='do not center the protein at the origin of the coordinate system')
    parser.add_argument('-ror',default='None,None',metavar='<atomselect string>,<atomselect string>',help='two comma-separated, single-quoted atomselect strings to define two groups of atoms whose centers of mass are aligned against the global z-axis')
    parser.add_argument('-v','--verbosity',action='count',help='output verbosity')

    args=parser.parse_args()
    
    if args.verbosity!=None:
        print('### Vebosity level: {}'.format(args.verbosity))
    else:
       args.verbosity=0
    Mut=MrgCmdLineAndFileContents(args.mut,args.mutfile,Mutation)
    Clv=MrgCmdLineAndFileContents(args.clv,args.clvfile,Cleavage)
    Gra=MrgCmdLineAndFileContents(args.gra,args.grafile,Graft)
    Att=MrgCmdLineAndFileContents(args.att,args.attfile,Attach)
    Uss=MrgCmdLineAndFileContents(args.ssbond,args.ssfile,SSBond)
    Usl=MrgCmdLineAndFileContents(args.link,args.linkfile,Link)
    UIC=args.ignore
    if len(args.topo)>0:
        CTopo.extend(args.topo)
    prefix=args.prefix
    PostMod['do_loop_mc']=args.rmi
    PostMod['do_gly_mc']=args.grel
    PostMod['Crot']=MrgCmdLineAndFileContents(args.crot,args.crotfile,Crot)
    PostMod['log_dcd_file']=args.logdcd
    PostMod['log_every']=args.logevery
    fixConflicts=not args.kc
    fixEngineeredMutations=args.rem
    psfgen=args.psfgen
    PostMod['center_protein']=~(args.noc)
    if args.ror!='None,None':
        PostMod['reorient_protein']=True
        PostMod['center_protein']=True
        PostMod['reorselstr']=args.ror.split(',')

    PDBfiles=args.pdbcif
    Molecules=[]
    if '.cif' in PDBfiles[0]:
        Molecules.append(Molecule(cif=PDBfiles[0],userLinks=Usl))
    elif '.pdb' in PDBfiles[0]:
        Molecules.append(Molecule(pdb=PDBfiles[0],userLinks=Usl))
    for p in PDBfiles[1:]:
        if '.cif' in p:
            Molecules.append(Molecule(cif=p))
        elif '.pdb' in p:
            Molecules.append(Molecule(pdb=p))
    Base=Molecules[0]
    Base.summarize()

    psfgen_fp=open(psfgen,'w')
    psfgen_fp.write('### This is an automatically generated psfgen input file\n')
    psfgen_fp.write('### created using cfapdbparse.py on {}\n'.format(date.today()))
    psfgen_fp.write('### cfapdbparse.py is part of the psfgen repository\n')
    psfgen_fp.write('### github.com:cameronabrams/psfgen/scripts\n')
    psfgen_fp.write('### questions to cfa22@drexel.edu\n')
    psfgen_fp.write('### command: python3 ')
    for a in sys.argv:
        psfgen_fp.write('{} '.format(a))
    psfgen_fp.write('\n')
    
    WriteHeaders(psfgen_fp,CTopo,LocTopo)
 
    if len(Clv)>0:
        Base.CleaveChains(Clv)

    Loops=Base.WritePsfgenInput(psfgen_fp,userMutations=Mut,fixConflicts=fixConflicts,fixEngineeredMutations=fixEngineeredMutations,prefix=prefix,userGrafts=Gra,userAttach=Att,userSSBonds=Uss,userIgnoreChains=UIC,removePDBs=True)

    ''' identify glycan segments '''
    glycan_segs=[]
    for c in Base.Chains.values():
        for s in c.Segments:
            if s.segtype=="GLYCAN":
                glycan_segs.append(s.segname)

    ''' fix crot replicas '''
    newcrots=[]
    for b in Base.Biomolecules:
        for t in b.biomt:
            if not t.isidentity():
                for c in PostMod['Crot']:
                    newcrots.append(c.replicate(t.get_replica_chainID(c.chainID)))
    PostMod['Crot'].extend(newcrots)
    post_pdb=WritePostMods(psfgen_fp,Base.psf_outfile,Base.pdb_outfile,PostMod,Loops,glycan_segs)

    Base.Tcl_PrependHeaderToPDB(post_pdb,psfgen_fp)

    psfgen_fp.write('exit\n')
    psfgen_fp.write('### thank you for using cfapdbparse.py!\n')
    print('"vmd -dispdev text -e {}" will generate {}/{}'.format(psfgen,Base.psf_outfile,post_pdb))
    
