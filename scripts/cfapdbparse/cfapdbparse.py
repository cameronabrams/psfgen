"""
    Parses PDB/mmCIF file(s) to build input file for VMD/psfgen
    Cameron F Abrams
    cfa22@drexel.edu

"""

import sys
import operator
import argparse
import os
import random
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
    lobsavevery=-1
    if 'log_dcd_file' in PostMod:
        logfile=PostMod['log_dcd_file']
    if 'log_every' in PostMod:
        logevery=PostMod['log_every']
    if 'log_save_every' in PostMod:
        logsaveevery=PostMod['log_save_every']
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
    if 'do_preheal_min_smd' in PostMod and PostMod['do_preheal_min_smd']:
        # 1. VMD: For each loop, lay down

        for l in sorted(Loops, key=lambda x: len(x.residues)):
            if (l.term and len(l.residues)>2):
                fp.write('lay_loop $molid {} [range {} {} 1] {}\n'.format(l.replica_chainID,l.residues[0].resseqnum,l.residues[-1].resseqnum,100))
        # 2. NAMD: Minimize
        # 3. VMD: Check for pierced rings
        # 4. NAMD: SMD loop closure
        #     - generate cv.inp
        # 5. VMD: Ligate loop-frag peptide bonds
        # 6. NAMD: Minimize (this is already the next step in the workflow)
        
    if 'do_multiflex_mc' in PostMod and PostMod['do_multiflex_mc']:
        nc=1000
        rcut=4.0
        sigma=1.8
        epsilon=0.5
        mctemperature=3.0
        mck=10.0
        dstop=2.0
        sstop=2.0
        maxanglestep=60.0 # degrees
        if 'multiflex_mc_params' in PostMod:
            p=PostMod['multiflex_mc_params']
            nc=nc if 'maxcycles' not in p else p['maxcycles']
            rcut=rcut if 'rcut' not in p else p['rcut']
            sigma=sigma if 'sigma' not in p else p['sigma']
            epsilon=epsilon if 'epsilon' not in p else p['epsilon']
            mctemperature=mctemperature if 'temperature' not in p else p['temperature']
            mck=mck if 'k' not in p else p['k']
            dstop=dstop if 'dstop' not in p else p['dstop']
            sstop=sstop if 'sstop' not in p else p['sstop']
            maxanglestep=maxanglestep if 'maxanglestep' not in p else p['maxanglestep']
        fp.write('set mcp [dict create]\n')
        fp.write('dict set mcp nc {}\n'.format(nc))
        fp.write('dict set mcp rcut {}\n'.format(rcut))
        fp.write('dict set mcp sigma {}\n'.format(sigma))
        fp.write('dict set mcp epsilon {}\n'.format(epsilon))
        fp.write('dict set mcp temperature {}\n'.format(mctemperature))
        fp.write('dict set mcp mck {}\n'.format(mck))
        fp.write('dict set mcp dstop {}\n'.format(dstop))
        fp.write('dict set mcp sstop {}\n'.format(sstop))
        fp.write('dict set mcp maxanglestep {}\n'.format(maxanglestep))
        fp.write('set bg [atomselect $molid "noh"]\n')
 #       fp.write('set loopindex 0\n')
 #       fp.write('set loops {\n')
        # build rotsel as as all atom indices in selection with rotatable bonds
        #  that is all atoms in all residues except for the C and O of last residue in each loop
        loopsel_substr=[]
        fa_substr=[]
        ca_substr=[]
        c_substr=[]
        #Loops.sort(key=lambda l: len(l.residues))
        for l in Loops:
            if l.term and len(l.residues)>1:
 #               fp.write('{{ {} {} {} }}\n'.format(l.replica_chainID,l.residues[0].resseqnum,l.residues[-1].resseqnum))
                loopsel_substr.append(' (chain {} and resid {} to {} and not (resid {} and name C O) )'.format(l.replica_chainID,l.residues[0].resseqnum,l.residues[-1].resseqnum,l.residues[-1].resseqnum))
                fa_substr.append(' (chain {} and resid {} and name CA) '.format(l.replica_chainID,l.residues[0].resseqnum))
                ca_substr.append(' (chain {} and resid {} and name CA) '.format(l.replica_chainID,l.residues[-1].resseqnum))
                c_substr.append(' (chain {} and resid {} and name C) '.format(l.replica_chainID,l.residues[-1].resseqnum))
        loopsel=' or '.join(loopsel_substr)
        fa_sel=' or '.join(fa_substr)
        ca_sel=' or '.join(ca_substr)
        c_sel=' or '.join(c_substr)
        loopsel='(protein and ('+loopsel+'))'
        fa_sel='(protein and ('+fa_sel+'))'
        ca_sel='(protein and ('+ca_sel+'))'
        c_sel='(protein and ('+c_sel+'))'
        fp.write('set fa [[atomselect $molid "{}"] get index]\n'.format(fa_sel))
        fp.write('set i [[atomselect $molid "{}"] get index]\n'.format(ca_sel))
        fp.write('set j [[atomselect $molid "{}"] get index]\n'.format(c_sel))

        if len(GlycanSegs)>0:
            glysel='(segname '+' '.join(GlycanSegs)+')'
            rotsel=loopsel+' or '+glysel
            fp.write('set gra {}\n')
            fp.write('set gi {}\n')
            fp.write('set gj {}\n')
            fp.write('set glycan_segs [list '+' '.join(GlycanSegs)+']\n')
            fp.write('set ng [llength $glycan_segs]\n')
            fp.write('foreach g $glycan_segs {\n')
            fp.write('   set sel [atomselect $molid "segname $g"]\n')
            fp.write('   set rid [$sel get resid]\n')
            fp.write('   set root [lindex [lsort -unique -real $rid] 0]\n')
            fp.write('   lappend gra [[atomselect $molid "segname $g and name C1 and resid $root"] get index]\n')
            fp.write('   lappend gi -1\n')
            fp.write('   lappend gj -1\n')
            fp.write('}\n')
            fp.write(r'set fa [list {*}$fa {*}$gra]'+'\n')
            fp.write(r'set i [list {*}$i {*}$gi]'+'\n')
            fp.write(r'set j [list {*}$j {*}$gj]'+'\n')

        fp.write('set rotsel [atomselect $molid "{}"]\n'.format(rotsel))
        fp.write('dict set atomind fa $fa\n'.format(fa_sel))
        fp.write('dict set atomind i $i\n'.format(ca_sel))
        fp.write('dict set atomind j $j\n'.format(c_sel))
        fp.write('do_multiflex_mc $molid $rotsel atomind mcp [irand_dom 1000 9999] $logid {} {}\n'.format(logevery,logsaveevery))
 #       fp.write('           }\n')
 #       fp.write('set nloops [llength $loops]\n')
        # build fa as list of all fixed atoms

        # build ca and c as lists of all ca-c indices (-1,-1 for glycans)

    
    new_pdb_out=prefix+'_mod.pdb'
    fp.write('$a writepdb {}\n'.format(new_pdb_out))
    if logdcd:
        fp.write('set loga [atomselect $logid all]\n')
        fp.write('animate write dcd {} waitfor all sel $loga $logid\n'.format(logfile))
        fp.write('mol delete $logid\n')
    return new_pdb_out

def WriteHeaders(fp,charmm_topologies,local_topologies):
    fp.write('#### BEGIN HEADER\n')
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
    fp.write('#### END HEADER\n')

def MrgCmdLineAndFileContents(cl_list,filename,typ):
    if filename!='':
        with open(filename,'r') as f:
           for l in f:
               if l[0]!='#':
                   cl_list.append(typ(l))
    return cl_list

def DictFromString(string):
    #print('parsing {}'.format(string))
    my_dict = {}
    if len(string)>0:
        items=string.split(',')
        for i in items:
            kv=i.split('=')
            k=kv[0]
            v=kv[1]
            my_dict[k]=v
    return my_dict

if __name__=='__main__':
    seed=random.randint(0,100000)
    temperature=400
    nummin=1000
    numsteps=2000
    target_numsteps=20000
    parser=argparse.ArgumentParser()
    print('cfapdbparse {} / python {}'.format(date.today(),sys.version.replace('\n',' ').split(' ')[0]))
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
    parser.add_argument('-logsaveevery',metavar='<int>',default=1,help='number of MC accepts between log writes to disk')
 #   parser.add_argument('-rlxloops',action='store_true',help='asks psfgen to use the loopMC module to relax modeled-in loops of residues missing from PDB')
    parser.add_argument('-rlxmc',action='store_true',help='asks psfgen to use do_multiflex_mc module to relax modeled-in loops of residues missing from PDB and glycans')
#    parser.add_argument('-loopmcparams',metavar='<param1=val1,param2=val2,...>',default='',help='Loop Monte Carlo parameters')
    parser.add_argument('-rlxmcparams',metavar='<param1=val1,param2=val2,...>',default='',help='Loop Monte Carlo parameters')
 #   parser.add_argument('-rlxgly',action='store_true',help='asks psfgen to use the loopMC module to relax modeled-in glycans missing from PDB')
 #   parser.add_argument('-glymcparams',metavar='<param1=val1,param2=val2,...>',default='',help='Glycan Monte Carlo parameters')
    parser.add_argument('-smdheal',action='store_true',help='asks psfgen to prep for a healing MD simulations to close missing loops')
    parser.add_argument('-kc',action='store_true',help='ignores SEQADV records indicating conflicts; if unset, residues in conflict are mutated to their proper identities')
    parser.add_argument('-rem',action='store_true',help='revert engineered mutations listed in SEQADV records')
    parser.add_argument('-noc',action='store_true',help='do not center the protein at the origin of the coordinate system')
    parser.add_argument('-ror',default='None,None',metavar='<atomselect string>,<atomselect string>',help='two comma-separated, single-quoted atomselect strings to define two groups of atoms whose centers of mass are aligned against the global z-axis')
    parser.add_argument('-v','--verbosity',action='count',help='output verbosity')
    parser.add_argument('-postscript',metavar='<name>',default='postscript.sh',help='autogenerated shell script to be run after')
    parser.add_argument('-pe',metavar='<int>',default=8,type=int,help='number of processors to indicated in NAMD inputs')

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
#    PostMod['do_loop_mc']=args.rlxloops
#    PostMod['loop_mc_params']=DictFromString(args.loopmcparams)
#    PostMod['do_gly_mc']=args.rlxgly
#    PostMod['gly_mc_params']=DictFromString(args.glymcparams)
    PostMod['do_multiflex_mc']=args.rlxmc
    PostMod['multiflex_mc_params']=DictFromString(args.rlxmcparams)
    PostMod['do_preheal_min_smd']=args.smdheal
    PostMod['Crot']=MrgCmdLineAndFileContents(args.crot,args.crotfile,Crot)
    PostMod['log_dcd_file']=args.logdcd
    PostMod['log_every']=args.logevery
    PostMod['log_save_every']=args.logsaveevery
    fixConflicts=not args.kc
    fixEngineeredMutations=args.rem
    psfgen=args.psfgen
    PostMod['center_protein']=~(args.noc)
    if args.ror!='None,None':
        PostMod['reorient_protein']=True
        PostMod['center_protein']=True
        PostMod['reorselstr']=args.ror.split(',')
    postscriptname=args.postscript
    npe=args.pe
    #print('-pe {:d}; NAMD will use {:d} processors.'.format(args.pe,npe))
    namdp='+p{:d}'.format(npe)
 
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

    ''' this will issue the final 'writepsf' and 'writepdb commands '''
    Loops=Base.WritePsfgenInput(psfgen_fp,userMutations=Mut,fixConflicts=fixConflicts,
                               fixEngineeredMutations=fixEngineeredMutations,prefix=prefix,
                               userGrafts=Gra,userAttach=Att,userSSBonds=Uss,userIgnoreChains=UIC,
                               removePDBs=True)

    ''' PostMods alter coordinates to ease minimization; psf is not modified further
        Regardless of whether any modifications are done or not, this will always write 
        a *_mod.pdb coordinate file '''

    ''' identify glycan segments '''
    glycan_segs=[]
    for c in Base.Chains.values():
        for s in c.Segments:
            if s.segtype=="GLYCAN":
                glycan_segs.append(s.segname)
    ''' generate crot replicas '''
    newcrots=[]
    for b in Base.Biomolecules:
        for t in b.biomt:
            if not t.isidentity():
                for c in PostMod['Crot']:
                    newcrots.append(c.replicate(t.get_replica_chainID(c.chainID)))
    PostMod['Crot'].extend(newcrots)
    post_pdb=WritePostMods(psfgen_fp,Base.psf_outfile,Base.pdb_outfile,PostMod,Loops,glycan_segs)
    # save useful header records to the output pdb
    Base.Tcl_PrependHeaderToPDB(post_pdb,psfgen_fp)

    psfgen_fp.write('exit\n')
    psfgen_fp.write('### thank you for using cfapdbparse.py!\n')

    ''' Generate the postscript '''
    print('Run the script {} to complete the build.'.format(postscriptname))
    print('After running {}, "read CURRPSF CURRPDB < .tmpvar" will set those variables.'.format(postscriptname))
    print('cfapdbpyparse ends.')
    fp=open(postscriptname,'w')
    fp.write(r'#!/bin/bash'+'\n')
    fp.write('# {}: completes the build of {}\n'.format(postscriptname,Base.psf_outfile))
    fp.write('TASK=$1\n')
    fp.write('echo "Completing the task-'+r'${TASK}'+' build of {}"\n'.format(Base.psf_outfile))
    fp.write('echo "VMD/PSFGEN) script={}'.format(psfgen))
    fp.write(' log=psfgen'+r'${TASK}'+'.log')
    fp.write(' psf={}'.format(Base.psf_outfile))
    fp.write(' pdb={}"\n'.format(post_pdb))
    fp.write(r'$VMD -dispdev text -e '+'{}'.format(psfgen)+r' 2&> psfgen${TASK}.log'+'\n')
    # save the patches!
    fp.write("cat {} | sed \'1,/#### BEGIN PATCHES/d;/#### END PATCHES/,$d\' > patches.inp\n".format(psfgen))
    fp.write('echo "structure {}" > tmpnamdheader\n'.format(Base.psf_outfile))
    fp.write('echo "coordinates {}" >> tmpnamdheader\n'.format(post_pdb))
    fp.write('cat tmpnamdheader $PSFGEN_BASEDIR/templates/vac.namd | sed s/%NUMMIN%/{}/ |'.format(nummin))
    fp.write(' sed s/%NUMSTEPS%/{}/ | sed s/%OUT%/tmpconfig/g | sed s/%SEED%/{}/g |'.format(numsteps,random.randint(0,10000)))
    fp.write(' sed s/%TEMPERATURE%/{}/g'.format(temperature))
    fp.write(r' > run${TASK}-1.namd'+'\n')
    fp.write('rm tmpnamdheader\n')
    fp.write('echo "NAMD2) {} config=run'.format(namdp)+r'${TASK}'+'-1.run log=run'+r'${TASK}'+'-1.log outputname=tmpconfig"\n')
    fp.write(r'$CHARMRUN '+namdp+r' $NAMD2 run${TASK}-1.namd > run${TASK}-1.log'+'\n')
    fp.write('if [ $? -ne 0 ]; then\n')
    fp.write('   echo "NAMD failed.  Check log file run'+r'${TASK}'+'-1.log."\n')
    fp.write('   exit\n')
    fp.write('fi\n')
    fp.write('echo "VMD) script=\$PSFGEN_BASEDIR/scripts/namdbin2pdb.tcl args={} tmpconfig.coor tmp.pdb'.format(Base.psf_outfile))
    fp.write(' log=namdbin2pdb'+r'${TASK}'+'-1.log"\n')
    fp.write(r'$VMD -dispdev text -e $PSFGEN_BASEDIR/scripts/namdbin2pdb.tcl -args '+'{} tmpconfig.coor tmp.pdb 2&> namdbin2pdb'.format(Base.psf_outfile)+r'${TASK}'+'-1.log\n')
    fp.write('cat charmm_header.pdb tmp.pdb > config.pdb\n')
    #fp.write('rm charmm_header.pdb tmp.pdb\n')
    fp.write('echo "VMD) script=$PSFGEN_BASEDIR/scripts/ringp.tcl args={} {}'.format(Base.psf_outfile,'config.pdb'))
    fp.write(' log=ringp'+r'${TASK}'+'-1.log"\n')
    fp.write(r'$VMD -dispdev text -e $PSFGEN_BASEDIR/scripts/ringp.tcl -args '+'{} {} 2&> ringp'+r'${TASK}'+'-1.log\n')
    fp.write('npiercings=`grep -c pierces ringp'+r'${TASK}'+'-1.log`\n')
    fp.write(r'if [[ $npiercings -gt 0 ]]; then'+'\n')
    fp.write(r'  echo "Error: There are $npiercings piercings in '+'{}"\n'.format('config.pdb'))
    fp.write('  grep pierces ringp'+r'${TASK}'+'-1.log\n')
    fp.write('  echo "Change your relaxation parameters and try again."\n')
    fp.write('  exit\n')
#    fp.write('else\n')
#    fp.write('  echo "No pierced rings found."\n')
    fp.write('fi\n')
    if 'do_preheal_min_smd' in PostMod and PostMod['do_preheal_min_smd']:
        fp.write('cat > heal_these.inp << EOF\n')
        for l in sorted(Loops, key=lambda x: len(x.residues)):
            if (l.term and len(l.residues)>2):
                #fp.write('# will try to heal bond between {} and {} on chain {}...\n'.format(l.residues[-1].resseqnum,l.nextfragntermres,l.replica_chainID))
                fp.write('{} {} {}\n'.format(l.replica_chainID,l.residues[-1].resseqnum,l.nextfragntermres))
#                fp.write('lay_loop $molid {} [range {} {} 1] {}\n'.format(l.replica_chainID,l.residues[0].resseqnum,l.residues[-1].resseqnum,100))
        fp.write('EOF\n')
        # measures to find the initial distances; generated fixed.pdb to fix the N atoms 
        fp.write('echo "VMD) script=\$PSFGEN_BASEDIR/scripts/measure_bonds.tcl args={} {} heal_these.inp'.format(Base.psf_outfile,'config.pdb'))
        fp.write(' log=heal'+r'${TASK}'+'.log"\n')
        fp.write(r'$VMD -dispdev text -e $PSFGEN_BASEDIR/scripts/measure_bonds.tcl -args ')
        fp.write('{} {} heal_these.inp 2&> heal'.format(Base.psf_outfile,'config.pdb')+r'${TASK}'+'.log\n')
        fp.write('if [ -f cv.inp ]; then rm cv.inp; fi\n')
        fp.write('touch cv.inp\n')
        fp.write('while IFS=" " read -r C L R B; do\n')
        fp.write(r'  cat $PSFGEN_BASEDIR/templates/cv-template.in | sed s/%C%/$C/g |')
        fp.write(r' sed s/%NAME%/${C}${L}/g | sed s/%I%/$L/g | sed s/%J%/$R/g | sed s/%R0%/$B/g |')
        fp.write(' sed s/%TARGETNUMSTEPS%/{}/ >> cv.inp ;\n'.format(target_numsteps))
        fp.write('done < heal_these.inp\n')
        fp.write('echo "structure {}" > tmpnamdheader\n'.format(Base.psf_outfile))
        fp.write('echo "coordinates {}" >> tmpnamdheader\n'.format('config.pdb'))
        fp.write('cat tmpnamdheader $PSFGEN_BASEDIR/templates/vac.namd |')
        fp.write(' sed s/%NUMMIN%/{}/ | sed s/%NUMSTEPS%/{}/ |'.format(0,int(1.5*target_numsteps)))
        fp.write(' sed s/%OUT%/tmpconfig/g | sed s/%SEED%/{}/g |'.format(random.randint(0,10000)))
        fp.write(' sed s/%TEMPERATURE%/{}/g |'.format(temperature))
        fp.write(r' sed "41 i fixedatoms on" |')
        fp.write(r' sed "42 i fixedatomsfile fixed.pdb" |')
        fp.write(r' sed "43 i fixedatomscol B" |')
        fp.write(r' sed "44 i colvars on" |')
        fp.write(r' sed "45 i colvarsconfig cv.inp" ')
        fp.write(' > run'+r'${TASK}'+'-2.namd\n')
        fp.write('rm tmpnamdheader\n')
        fp.write('echo "NAMD2) config=run'+r'${TASK}'+'-2.namd log=run'+r'${TASK}'+'-2.log outputname=tmpconfig"\n')
        fp.write(r'$CHARMRUN '+namdp+' $NAMD2 run'+r'${TASK}'+'-2.namd > run${TASK}-2.log'+'\n')
        fp.write('if [ $? -ne 0 ]; then\n')
        fp.write('   echo "NAMD failed.  Check log file run'+r'${TASK}'+'-2.log."\n')
        fp.write('   exit\n')
        fp.write('fi\n')
        fp.write('echo "VMD) script=\$PSFGEN_BASEDIR/scripts/namdbin2pdb.tcl args={} tmpconfig.coor tmpconfig2.pdb'.format(Base.psf_outfile))
        fp.write(' log=namdbin2pdb'+r'${TASK}'+'-2.log"\n')
        fp.write(r'$VMD -dispdev text -e $PSFGEN_BASEDIR/scripts/namdbin2pdb.tcl -args '+'{} tmpconfig.coor tmpconfig2.pdb 2&> namdbin2pdb'.format(Base.psf_outfile)+r'${TASK}'+'-2.log\n')
        # prepend the charmm header to the pdb file
        fp.write('cat charmm_header.pdb tmp.pdb > config2.pdb\n')
        fp.write('cat > the_healing_patches.inp << EOF\n')
        for l in sorted(Loops, key=lambda x: len(x.residues)):
            if (l.term and len(l.residues)>2):
                #fp.write('# will try to heal bond between {} and {} on chain {}...\n'.format(l.residues[-1].resseqnum,l.nextfragntermres,l.replica_chainID))
                fp.write('patch HEAL {c}:{ll} {c}:{l} {c}:{r} {c}:{rr}\n'.format(c=l.replica_chainID,
                            ll=l.residues[-2].resseqnum,l=l.residues[-1].resseqnum,r=l.nextfragntermres,rr=(l.nextfragntermres+1)))
        fp.write('EOF\n')
        fp.write('cat $PSFGEN_BASEDIR/scripts/ligations.tcl | sed "/#### LIGATION LIST STARTS/r the_healing_patches.inp"  > do_the_healing.tcl\n')
        fp.write('echo "VMD/PSFGEN) script=do_the_healing.tcl args={} {} {} {}'.format(Base.psf_outfile,
        'tmpconfig2.pdb','ligated.psf','tmp2config2.pdb'))
        fp.write(' log=ligations'+r'${TASK}'+'.log"\n')
        fp.write(r'$VMD -dispdev text -e do_the_healing.tcl -args '+'{} {} {} {} 2&> ligations'.format(Base.psf_outfile,
        'tmpconfig2.pdb','ligated.psf','tmp2config2.pdb')+r'${TASK}'+'.log\n')
        fp.write('echo "structure {}" > tmpnamdheader\n'.format('ligated.psf'))
        fp.write('echo "coordinates {}" >> tmpnamdheader\n'.format('tmp2config2.pdb'))
        fp.write('cat tmpnamdheader $PSFGEN_BASEDIR/templates/vac.namd |')
        fp.write(' sed s/%NUMMIN%/{}/ | sed s/%NUMSTEPS%/{}/ |'.format(nummin,numsteps))
        fp.write(' sed s/%OUT%/tmpconfig3/g | sed s/%SEED%/{}/g |'.format(random.randint(0,10000)))
        fp.write(' sed s/%TEMPERATURE%/{}/g '.format(temperature))
        fp.write(' > run'+r'${TASK}'+'-3.namd\n')
        fp.write('rm tmpnamdheader\n')
       # fp.write('echo "Running namd2 min on ligated system {} {}; output in run'.format('ligated.psf','tmp2config2.pdb')+r'${TASK}'+'-3.log"\n')
        fp.write('echo "NAMD2) config=run'+r'${TASK}'+'-3.namd log=run'+r'${TASK}'+'-3.log outputname=tmpconfig3"\n')
        fp.write(r'$CHARMRUN '+namdp+' $NAMD2 run'+r'${TASK}'+'-3.namd > run'+r'${TASK}'+'-3.log'+'\n')
        fp.write('if [ $? -ne 0 ]; then\n')
        fp.write('   echo "NAMD failed.  Check log file run'+r'${TASK}'+'-3.log."\n')
        fp.write('   exit\n')
        fp.write('fi\n')
        fp.write('echo "VMD) script=\$PSFGEN_BASEDIR/scripts/namdbin2pdb.tcl args={} tmpconfig3.coor tmpconfig4.pdb'.format('ligated.psf'))
        fp.write(' log=namdbin2pdb'+r'${TASK}'+'-3.log"\n')
        fp.write(r'$VMD -dispdev text -e $PSFGEN_BASEDIR/scripts/namdbin2pdb.tcl -args '+'{} tmpconfig3.coor tmpconfig4.pdb 2&> namdbin2pdb'.format('ligated.psf')+r'${TASK}'+'-3.log\n')
        fp.write('cat charmm_header.pdb tmpconfig4.pdb > config2.pdb\n')
        fp.write('echo "VMD) script=$PSFGEN_BASEDIR/scripts/ringp.tcl args={} {}'.format('ligated.psf','config2.pdb'))
        fp.write(' log=ringp'+r'${TASK}'+'-2.log"\n')
        #fp.write('echo "Checking {} for pierced rings; log is ringp'.format('config2.pdb')+r'${TASK}'+'-2.log"\n')
        fp.write(r'$VMD -dispdev text -e $PSFGEN_BASEDIR/scripts/ringp.tcl -args '+'{} {} 2&> ringp'+r'${TASK}'+'-2.log\n'.format('ligated.psf','config2.pdb'))
        fp.write('npiercings=`grep -c pierces ringp'+r'${TASK}'+'-2.log`\n')
        fp.write(r'if [[ $npiercings -gt 0 ]]; then'+'\n')
        fp.write(r'  echo "Error: There are $npiercings piercings in '+'{}"\n'.format('config2.pdb'))
        fp.write('  grep pierces ringp'+r'${TASK}'+'-2.log\n')
        fp.write('  echo "Change your relaxation parameters and try again."\n')
        fp.write('  exit\n')
#        fp.write('else\n')
#        fp.write('  echo "No pierced rings found."\n')
        fp.write('fi\n')
        fp.write('echo {} {} > .tmpvar\n'.format('ligated.psf','config2.pdb'))
    else:
        fp.write('echo {} {} > .tmpvar\n'.format(Base.psf_outfile,'config2.pdb'))
    fp.write('# {} finishes.\n'.format(postscriptname))
    fp.close()
    os.system('chmod 744 {}'.format(postscriptname))

    
