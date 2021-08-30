"""
    Parses PDB/mmCIF file(s) to build input Tcl file for VMD/psfgen
    Cameron F Abrams
    cfa22@drexel.edu

"""

import sys
import operator
import argparse
import os
import random
import math
from datetime import date, datetime
from molecule import Molecule
from cleavage import Cleavage
from mutation import Mutation
from ssbond import SSBond
from graft import Graft
from crot import Crot
from attach import Attach
from link import Link
from deletion import Deletion
from modsfile import ModsFile
from atom import _PDBAtomNameDict_
from residue import Residue, _PDBResName123_, _pdb_glycans_, _pdb_ions_, _ResNameDict_PDB_to_CHARMM_, _ResNameDict_CHARMM_to_PDB_, get_residue

def vmd_instructions(fp,script,logname='tmp.log',args='',msg=''):
    fp.write('echo "VMD) script={} log={} msg: {}"\n'.format(script,logname,msg))
    if args!='':
        fp.write(r'$VMD -dispdev text -e '+script+r' -args '+args+r' > '+logname+' 2>&1\n')
    else:
        fp.write(r'$VMD -dispdev text -e '+script+r' > '+logname+' 2>&1\n')
    fp.write('if [ $? -ne 0 ]; then\n')
    fp.write('   echo "VMD failed.  Check the log file {}. Exiting."\n'.format(logname))
    fp.write('   exit 1\n')
    fp.write('fi\n')

def namd_instructions(fp,cfgname,psf,coor,outname,logname,
                      npe=8,numminsteps=0,numsteps=0,seed=0,template='vac.namd',
                      temperature=310,extras=[],msg='',stdparamfiles=[],localparamfiles=[],
                      stdcharmmdir=r'\$env(HOME)/charmm/toppar',
                      localcharmmdir=r'\$env(PSFGEN_BASEDIR)/charmm'):
    fp.write('cat $PSFGEN_BASEDIR/templates/{}'.format(template))
    fp.write('  | sed s/%OUT%/{}/g'.format(outname))
    fp.write('  | sed s/%NUMMIN%/{}/'.format(numminsteps))
    fp.write('  | sed s/%NUMSTEPS%/{}/'.format(numsteps))
    fp.write('  | sed s/%SEED%/{}/g'.format(seed))
    fp.write('  | sed s/%TEMPERATURE%/{}/g'.format(temperature))
    fp.write('  | sed "/#### SYSTEM CONFIGURATION FILES END/i structure {}"'.format(psf))
    fp.write('  | sed "/#### SYSTEM CONFIGURATION FILES END/i coordinates {}"'.format(coor))
    sentinelline='#### PARAMETER FILES END'
    for st in stdparamfiles:
        fp.write(' | sed "/{}/i parameters {}/{}" '.format(sentinelline,stdcharmmdir,st))
    for st in localparamfiles:
        fp.write(' | sed "/{}/i parameters {}/{}" '.format(sentinelline,localcharmmdir,st))
    sentinelline='#### EXTRAS END'
    for ex in extras:
        fp.write('  | sed "/'+sentinelline+'/i '+ex+'" ')
    fp.write(' > {}\n'.format(cfgname))
    namdp='+p{:d}'.format(npe)
    fp.write('echo "NAMD2) config={} log={} outputname={} msg={}"\n'.format(cfgname,logname,outname,msg))
    fp.write(r'$CHARMRUN '+namdp+r' $NAMD2 '+cfgname+r' > '+logname+'\n')
    fp.write('if [ $? -ne 0 ]; then\n')
    fp.write('   echo "NAMD failed.  Check log file {}. Exiting."\n'.format(logname))
    fp.write('   exit 1\n')
    fp.write('fi\n')

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
    if 'lay_cycles' in PostMod:
        lay_cycles=PostMod['lay_cycles']

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
            fp.write('set ca [measure center [atomselect $molid "protein and {}"] weight mass]\n'.format(PostMod['reorselstr'][0]))
            fp.write('set cb [measure center [atomselect $molid "protein and {}"] weight mass]\n'.format(PostMod['reorselstr'][1]))
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
        fp.write(crot.psfgen_str(molid=r'$molid'))
        if logdcd:
            fp.write('log_addframe $molid $logid\n')
    if 'do_preclose_min_smd' in PostMod and PostMod['do_preclose_min_smd']:
        lay_cycles=100
        if 'preclose_params' in PostMod:
            p=PostMod['preclose_params']
            lay_cycles=lay_cycles if 'lay_cycles' not in p else p['lay_cycles']
        for l in sorted(Loops, key=lambda x: len(x.residues)):
            if (l.term and len(l.residues)>2):
                fp.write('lay_loop $molid {} [range {} {} 1] {}\n'.format(l.replica_chainID,
                          l.residues[0].resseqnum,l.residues[-1].resseqnum,lay_cycles))

        
    if 'do_multiflex_mc' in PostMod and PostMod['do_multiflex_mc']:
        nc=1000
        rcut=4.0
        sigma=1.8
        epsilon=0.5
        cutoff=math.pow(2,(1./6.))*sigma
        shift=epsilon
        mctemperature=3.0
        mck=10.0
        dstop=2.0
        sstop=2.0
        maxanglestep=60.0 # degrees
        do_loops = 0
        do_gly = 1
        if 'multiflex_mc_params' in PostMod:
            p=PostMod['multiflex_mc_params']
            nc=nc if 'maxcycles' not in p else p['maxcycles']
            rcut=rcut if 'rcut' not in p else p['rcut']
            sigma=sigma if 'sigma' not in p else p['sigma']
            epsilon=epsilon if 'epsilon' not in p else p['epsilon']
            shift=shift if 'shift' not in p else p['shift']
            cutoff=cutoff if 'cutoff' not in p else p['cutoff']
            mctemperature=mctemperature if 'temperature' not in p else p['temperature']
            mck=mck if 'k' not in p else p['k']
            dstop=dstop if 'dstop' not in p else p['dstop']
            sstop=sstop if 'sstop' not in p else p['sstop']
            maxanglestep=maxanglestep if 'maxanglestep' not in p else p['maxanglestep']
            do_loops=do_loops if 'loops' not in p else p['loops']
            do_gly=do_gly if 'gly' not in p else p['gly']
        fp.write('set mcp [dict create]\n')
        fp.write('dict set mcp nc {}\n'.format(nc))
        fp.write('dict set mcp cellsize {}\n'.format(rcut))
        fp.write('dict set mcp ljsigma {}\n'.format(sigma))
        fp.write('dict set mcp ljepsilon {}\n'.format(epsilon))
        fp.write('dict set mcp ljcutoff {}\n'.format(cutoff))
        fp.write('dict set mcp ljshift {}\n'.format(shift))
        fp.write('dict set mcp temperature {}\n'.format(mctemperature))
        fp.write('dict set mcp mck {}\n'.format(mck))
        fp.write('dict set mcp dstop {}\n'.format(dstop))
        fp.write('dict set mcp sstop {}\n'.format(sstop))
        fp.write('dict set mcp maxanglestep {}\n'.format(maxanglestep))
        fp.write('set bg [atomselect $molid "noh"]\n')
        fp.write('set fa {}\n')
        fp.write('set i {}\n')
        fp.write('set j {}\n')
 #       fp.write('set loopindex 0\n')
 #       fp.write('set loops {\n')
        # build rotsel as as all atom indices in selection with rotatable bonds
        #  that is all atoms in all residues except for the C and O of last residue in each loop
        rotsel=''
        if len(Loops)>0 and do_loops == 1:
            loopsel_substr=[]
            fa_substr=[]
            ca_substr=[]
            c_substr=[]
            #Loops.sort(key=lambda l: len(l.residues))
            for l in Loops:
                if l.term and len(l.residues)>1:
 #                   fp.write('{{ {} {} {} }}\n'.format(l.replica_chainID,l.residues[0].resseqnum,l.residues[-1].resseqnum))
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
            fp.write('set lfa [[atomselect $molid "{}"] get index]\n'.format(fa_sel))
            fp.write('set lca [[atomselect $molid "{}"] get index]\n'.format(ca_sel))
            fp.write('set lc [[atomselect $molid "{}"] get index]\n'.format(c_sel))
            fp.write(r'set fa [list {*}$fa {*}$lfa]'+'\n')
            fp.write(r'set i [list {*}$i {*}$lca]'+'\n')
            fp.write(r'set j [list {*}$j {*}$lc]'+'\n')
            rotsel=loopsel
            
        if len(GlycanSegs)>0 and do_gly == 1:
            glysel='(segname '+' '.join(GlycanSegs)+')'
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

        if len(rotsel)>0:
            rotsel = rotsel+' or '+glysel
        else:
            rotsel=glysel

        fp.write('set rotsel [atomselect $molid "{}"]\n'.format(rotsel))
        fp.write('dict set atomind fa $fa\n')
        fp.write('dict set atomind i $i\n')
        fp.write('dict set atomind j $j\n')
        fp.write('do_multiflex_mc $molid $rotsel atomind mcp [irand_dom 1000 9999] $logid {} {}\n'.format(logevery,logsaveevery))
    
    new_pdb_out=prefix+'_mod.pdb'
    fp.write('$a writepdb {}\n'.format(new_pdb_out))
    if logdcd:
        fp.write('set loga [atomselect $logid all]\n')
        fp.write('animate write dcd {} waitfor all sel $loga $logid\n'.format(logfile))
        fp.write('mol delete $logid\n')
    return new_pdb_out

def CommonPSFGENheader(fp,charmm_topologies,local_topologies):
    for t in charmm_topologies:
        fp.write('topology $TOPPARDIR/{}\n'.format(t))
    for t in local_topologies:
        fp.write('topology $LOCAL_TOPPARDIR/{}\n'.format(t))

def WriteHeaders(fp,charmm_topologies,local_topologies,pdbaliases):
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
    fp.write('psfcontext mixedcase\n')
    CommonPSFGENheader(fp,charmm_topologies,local_topologies)
    #print(pdbaliases)
    for al in pdbaliases:
        fp.write('pdbalias {}\n'.format(al))
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

def DefOrDict(d,varname,default):
    return default if varname not in d else d[varname]

def GetStreamFileNames(giventopos):
    streamfilesnames=[]
    for t in giventopos:
        if t[-3:] == 'str':
            streamfilesnames.append(t)
    return streamfilesnames

if __name__=='__main__':
    seed=random.randint(0,100000)

    parser=argparse.ArgumentParser()
    print('cfapdbparse.py {} / python {}'.format(date.today(),sys.version.replace('\n',' ').split(' ')[0]))
    i=1
    Molecules=[]
    Mut=[]
    Clv=[]
    Uss=[]
    UIC=[]
    # defaults
    psfgen='mkpsf.tcl'
    CTopo=['top_all36_prot.rtf','top_all35_ethers.rtf','top_all36_cgenff.rtf','top_all36_lipid.rtf',
           'top_all36_na.rtf','stream/carb/toppar_all36_carb_glycopeptide.str']
    # default local topologies: these are specially modified charmm str files that get rid of things that PSFGEN can't handle
    LocTopo=['top_all36_carb.rtf','toppar_water_ions.str']
    StdParamFiles=['par_all36_prot.prm','par_all36_carb.prm','par_all36_lipid.prm','par_all36_na.prm','par_all36_cgenff.prm','stream/carb/toppar_all36_carb_glycopeptide.str']
    LocalParamFiles=['toppar_water_ions.str']
    PDBAliases=['residue HIS HSD','atom ILE CD1 CD','residue NAG BGNA','atom BGNA C7 C',
                        'atom BGNA O7 O','atom BGNA C8 CT','atom BGNA N2 N','residue SIA ANE5',
                        'atom ANE5 C10 C','atom ANE5 C11 CT','atom ANE5 N5 N','atom ANE5 O1A O11',
                        'atom ANE5 O1B O12','atom ANE5 O10 O','atom VCG C01 C1','atom VCG C01 C1','atom VCG C02 C2',
                        'atom VCG C03 C3','atom VCG C04 C4','atom VCG C05 C5','atom VCG C06 C6','atom VCG C07 C7',
                        'atom VCG C08 C8','atom VCG C09 C9','residue EIC LIN']
    PostMod={}
    PostMod['center_protein']=True
    prefix='x01_'
    fixConflicts=True
    PostMod['do_loop_mc']=False
    PostMod['do_gly_mc']=False
    PostMod['Crot']=[]

    parser.add_argument('pdbcif',nargs='+',metavar='<?.pdb|cif>',type=str,help='Name(s) of pdb or CIF file to parse; First is treated as the base molecule')
    parser.add_argument('-ba','--biological-assembly',metavar='#',default=0,type=int,help='Biological assembly to construct; one may be selected from those defined in PDB/mmCIF metadata.  If not specified, the explicit model is built.')
    parser.add_argument('-charmmtopo',metavar='<name> ...',nargs='+',default=[],help='Additional (standard) CHARMM topology files in your CHARMM directory')
    parser.add_argument('-loctopo',metavar='<name> ...',nargs='+',default=[],help='Additional (local) CHARMM topology files in the $PSFGEN_BASEDIR/charmm directory')
    parser.add_argument('-charmmparam',metavar='<name> ...',nargs='+',default=[],help='Additional (standard) CHARMM parameter files in your CHARMM directory')
    parser.add_argument('-locparam',metavar='<name> ...',nargs='+',default=[],help='Additional (local) CHARMM parameter files in the $PSFGEN_BASEDIR/charmm directory')
    parser.add_argument('-prefix',metavar='<str>',default='x01_',help='Output PDB/PSF prefix; each file name will have the format <prefix><pdbcode>.pdb/psf, where <pdbcode> is the 4-letter PDB code of the base molecule.')
    parser.add_argument('-psfgen',metavar='<name>',default='mkpsf.tcl',help='name of TcL script generated as input to VMD/psfgen')
    parser.add_argument('-ignore',metavar='X ...',nargs='+',default=[],type=str,help='Specify chain(s) to ignore.')
    parser.add_argument('-modsfile',metavar='<name> ...',nargs='+',default=[],type=ModsFile,help='One (or more) modifications file(s) to rule them all.')
    parser.add_argument('-mut',metavar='C_OrrrN [C_OrrrN] ...',nargs='+',default=[],type=Mutation,help='One or more point-mutation specifications.  Format: C is chainID, O is one-letter residue code to mutate FROM, rrr is sequence number (can be any number of digits), and N is one-letter residue code to mutate TO.  Multiple mutation instances can be specified with one -mut.  Mutations are automatically replicated if there are BIOMT transformations.')
    parser.add_argument('-mutfile',metavar='<name>',default='',help='Input file listing mutation specifications')
    parser.add_argument('-delete',metavar='C_Orrr [C_Orrr] ...',nargs='+',default=[],type=Deletion,help='One or more single-residue deletion specifications.  Format: C is chainID, O is one-letter residue code, rrr is sequence number (can be any number of digits).  Multiple deletion instances can be specified with one -del.  Deletions are automatically replicated if there are BIOMT transformations.')
    parser.add_argument('-deletefile',metavar='<name>',default='',help='Input file listing deletion specifications')
    parser.add_argument('-clv',metavar='PrrrC [PrrrC] ...',nargs='+',default=[],type=Cleavage,help='One or more cleavage-site specifications.  Format: P is parent chain ID, rrr is residue number immediately N-terminal to the cleavage site, and C is the daughter chain ID that will begin immediately C-terminal to cleavage site. Multiple cleavage instances can be specified after one -clv.')
    parser.add_argument('-clvfile',metavar='<name>',default='',help='input file listing all cleavages (as an alternative to issuing multiple -clv arguments)')
    parser.add_argument('-gra',metavar='<str>,A:XXX-YYY,ZZZ,C:BBB ...',nargs='+',default=[],type=Graft,help='One or more graft specifications; graft resids XXX-YYY of chain A in pdb <str> to chain C of base molecule by overlapping resid ZZZ of chain A of graft and resid BBB of chain C of base.  Grafts are automatically replicated if there are BIOMT transformations.')
    parser.add_argument('-grafile',metavar='<name>',default='',help='Input file listing all grafts (as an alternative to issuing multiple -gra arguments)')
    parser.add_argument('-att',metavar='<str>,A:XXX-YYY,ZZZ,B:QQQ,C:BBB ...',nargs='+',default=[],type=Attach,help='One or more attachment specifications. Format: attach resids XXX-YYY of chain A using resid ZZZ (between XXX and YYY) in pdb <str> to chain C of base molecule at resid BBB by aligning resid QQQ of chain B from source to resid BBB of chain C of base.  CURRENTLY UNIMPLEMENTED!!!')
    parser.add_argument('-attfile',metavar='<name>',default='',help='Input file listing all attachments (as an alternative to issuing multiple -att arguments)')
    parser.add_argument('-crot',metavar='<str>,A,XXX[,YYY],### ...',default=[],nargs='+',type=Crot,help='One or more torsion rotation specifications. Format: <str> is one of phi, psi, omega, chi1, or chi2.  A is the chainID, XXX is the resid of owner of torson, and YYY (if given) marks the end of the sequence C-terminal to XXX that is reoriented by a backbone rotation. ### is the degrees of rotation.  C-rotations are automatically replicated if there are BIOMT transformations.')
    parser.add_argument('-crotfile',metavar='<name>',default='',help='Input file listing all torsion rotations requested (as an alternative to issuing multiple -crot arguments)')
    parser.add_argument('-ssbond',metavar='X_###-Y_### ...',default=[],nargs='+',type=SSBond,help='One or more disulfide bond specification(s) not already in input PDB/CIF: Format: X,Y are chainIDs and ### are resids; if residues are not CYS in wt or by mutations, there is no effect.  Because SSBonds can join chains together, they are NOT automatially replicated if there are BIOMT transformations.')
    parser.add_argument('-ssfile',metavar='<name>',default='',help='input file listing all disulfide bonds to add that are not already in the PDB file (as an alternative to issuing multiple -ssbond arguments)')
    parser.add_argument('-link',metavar='string',default=[],action='append',type=Link,help='PDB-format LINK record; must have exact spacing; multiple "-link" options can be supplied.')
    parser.add_argument('-linkfile',metavar='<name>',default='',help='Input file with PDB-format LINK records the user would like to enforce that are not in the PDB/CIF file')
    parser.add_argument('-pdbalias',metavar='<str>',default=[],nargs='+',help='One or more psfgen-formatted pdbalias with commas for spaces')
    parser.add_argument('-pdbaliasfile',metavar='<str>',default='',help='Input file containing psfgen-formatted pdbaliases')
    parser.add_argument('-logdcd',metavar='<name>.dcd',default='',help='Name of dcd logging file')
    parser.add_argument('-logevery',metavar='<int>',default=1,help='Number of MC accepts between successive frame logging')
    parser.add_argument('-logsaveevery',metavar='<int>',default=1,help='Number of MC accepts between log writes to disk')
 #   parser.add_argument('-rlxloops',action='store_true',help='asks psfgen to use the loopMC module to relax modeled-in loops of residues missing from PDB')
#    parser.add_argument('-rlxmc',action='store_true',help='Asks psfgen to use do_multiflex_mc module to relax modeled-in loops of residues missing from PDB and glycans')
#    parser.add_argument('-loopmcparams',metavar='<param1=val1,param2=val2,...>',default='',help='Loop Monte Carlo parameters')
#    parser.add_argument('-rlxmcparams',metavar='<param1=val1,param2=val2,...>',default='',help='Loop Monte Carlo parameters')
#    parser.add_argument('-rlxgly',action='store_true',help='asks psfgen to use the loopMC module to relax modeled-in glycans missing from PDB')
#    parser.add_argument('-glymcparams',metavar='<param1=val1,param2=val2,...>',default='',help='Glycan Monte Carlo parameters')
    parser.add_argument('-smdclose',action='store_true',help='Asks psfgen to prep for steered MD simulation to close loops')
    parser.add_argument('-smdcloseparams',metavar='<param1=val1,param2=val2,...>',default='',help='Parameters for steered MD for closing loops')
    parser.add_argument('-namdparams',metavar='<param1=val1,param2=val2,...>',default='',help='Parameters for NAMD runs')
    parser.add_argument('-kc',action='store_true',help='Ignores SEQADV records indicating conflicts between construct and sequence; if unset, residues in CONFLICT are mutated to their proper identities')
    parser.add_argument('-rem',action='store_true',default=False,help='Revert engineered mutations listed in SEQADV records')
    parser.add_argument('-noc',action='store_true',help='do not center the protein at the origin of the coordinate system')
    parser.add_argument('-ror',default='None,None',metavar='<atomselect string>,<atomselect string>',help='two comma-separated, single-quoted atomselect strings to define two groups of atoms whose centers of mass are aligned against the global z-axis')
    parser.add_argument('-v','--verbosity',action='count',help='output verbosity')
    parser.add_argument('-postscript',metavar='<name>',default='postscript.sh',help='autogenerated shell script to be run after vmd')
    parser.add_argument('-pe',metavar='<int>',default=8,type=int,help='number of processors to indicated in NAMD inputs')

    args=parser.parse_args()
    
    if len(sys.argv)>1:
        print('Command-line arguments: '+' '.join(sys.argv[1:]))

    if args.verbosity!=None:
        print('### Vebosity level: {}'.format(args.verbosity))
    else:
       args.verbosity=0

    if args.verbosity>0:
        for k,v in vars(args).items():
            if type(v) is list and len(v)>0:
                print('-{:s} '.format(k)+' '.join(v),end=' ')
            elif type(v) is str and len(v)>0:
                print('-{:s} {}'.format(k,v),end=' ')
            elif not type(v) is str and not type(v) is list:
                print('-{:s} {}'.format(k,v),end=' ')
        print()

    Mut=MrgCmdLineAndFileContents(args.mut,args.mutfile,Mutation)
    Clv=MrgCmdLineAndFileContents(args.clv,args.clvfile,Cleavage)
    Gra=MrgCmdLineAndFileContents(args.gra,args.grafile,Graft)
    Att=MrgCmdLineAndFileContents(args.att,args.attfile,Attach)
    Uss=MrgCmdLineAndFileContents(args.ssbond,args.ssfile,SSBond)
    Usl=MrgCmdLineAndFileContents(args.link,args.linkfile,Link)
    Del=MrgCmdLineAndFileContents(args.delete,args.deletefile,Deletion)
    if len(args.modsfile)>0:
        for mf in args.modsfile:
            if args.verbosity>0:
                mf.report()
            Mut.extend(mf.show_type(Mutation))
            Clv.extend(mf.show_type(Cleavage))
            Gra.extend(mf.show_type(Graft))
            Att.extend(mf.show_type(Attach))
            Uss.extend(mf.show_type(SSBond))
            Usl.extend(mf.show_type(Link))
            Del.extend(mf.show_type(Deletion))

    UPDBAliases=MrgCmdLineAndFileContents([' '.join(_.split(',')) for _ in args.pdbalias],args.pdbaliasfile,str)
    PDBAliases.extend(UPDBAliases)
    UIC=args.ignore
    CTopo.extend(args.charmmtopo)
    LocTopo.extend(args.loctopo)
    StdParamFiles.extend(args.charmmparam)
    LocalParamFiles.extend(args.locparam)

    prefix=args.prefix
#    PostMod['do_loop_mc']=args.rlxloops
#    PostMod['loop_mc_params']=DictFromString(args.loopmcparams)
 #   PostMod['do_gly_mc']=args.rlxgly
 #   PostMod['gly_mc_params']=DictFromString(args.glymcparams)
    #PostMod['do_multiflex_mc']=args.rlxmc
    #ostMod['multiflex_mc_params']=DictFromString(args.rlxmcparams)
    PostMod['do_preclose_min_smd']=args.smdclose
    PostMod['preclose_params']=DictFromString(args.smdcloseparams)
    PostMod['NAMD_params']=DictFromString(args.namdparams)
    PostMod['Crot']=MrgCmdLineAndFileContents(args.crot,args.crotfile,Crot)
    if len(args.modsfile)>0:
        for mf in args.modsfile:
            PostMod['Crot'].extend(mf.show_type(Crot))
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
    for k,v in PostMod.items():
        print('{}:'.format(k),v)

    postscriptname=args.postscript
    npe=args.pe
    #print('-pe {:d}; NAMD will use {:d} processors.'.format(args.pe,npe))
 
    PDBfiles=args.pdbcif
    Molecules=[]
    if '.cif' in PDBfiles[0]:
        Molecules.append(Molecule(cif=PDBfiles[0],userLinks=Usl,requestedBiologicalAssembly=args.biological_assembly))
    elif '.pdb' in PDBfiles[0]:
        Molecules.append(Molecule(pdb=PDBfiles[0],userLinks=Usl,requestedBiologicalAssembly=args.biological_assembly))
    for p in PDBfiles[1:]:
        if '.cif' in p:
            Molecules.append(Molecule(cif=p))
        elif '.pdb' in p:
            Molecules.append(Molecule(pdb=p))
    Base=Molecules[0]
    Base.summarize()

    psfgen_fp=open(psfgen,'w')
    psfgen_fp.write('### This is an automatically generated psfgen input file\n')
    psfgen_fp.write('### created using cfapdbparse.py on {} at {}\n'.format(date.today(),datetime.now().strftime('%H:%M:%S')))
    psfgen_fp.write('### cfapdbparse.py is part of the psfgen repository\n')
    psfgen_fp.write('### github.com:cameronabrams/psfgen/scripts\n')
    psfgen_fp.write('### questions to cfa22@drexel.edu\n')
    psfgen_fp.write('### command: python3 ')
    for a in sys.argv:
        psfgen_fp.write('{} '.format(a))
    psfgen_fp.write('\n')
    
    WriteHeaders(psfgen_fp,CTopo,LocTopo,PDBAliases)
 
    if len(Clv)>0:
        Base.CleaveChains(Clv)

    ''' this will issue the final 'writepsf' and 'writepdb commands '''
    Loops=Base.WritePsfgenInput(psfgen_fp,userMutations=Mut,userDeletions=Del,fixConflicts=fixConflicts,
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
    nummin=1000
    numsteps=2000
    temperature=310
    if 'NAMD_params' in PostMod:
        p=PostMod['NAMD_params']
        nummin=DefOrDict(p,'nummin',nummin)
        numsteps=DefOrDict(p,'numsteps',numsteps)
        temperature=DefOrDict(p,'temperature',temperature)
    currpdb=post_pdb
    currpsf=Base.psf_outfile
    print('Run the script {} to complete the build.'.format(postscriptname))
    print('After running {}, "read CURRPSF CURRPDB CURRCFG < .tmpvar" will set those variables.'.format(postscriptname))
    print('cfapdbpyparse ends.')
    fp=open(postscriptname,'w')
    fp.write(r'#!/bin/bash'+'\n')
    fp.write('# {}: completes the build of {}\n'.format(postscriptname,currpsf))
    fp.write('''
source $PSFGEN_BASEDIR/scripts/utils.sh
nesting_level=1
TASK=1
i=1
ARGC=$#
while [ $i -le $ARGC ] ; do
  if [ "${!i}" = "-task" ]; then
    i=$((i+1))
    TASK=${!i}
  elif [ "${!i}" = "-nesting-level" ]; then
    i=$((i+1))
    nesting_level=${!i}
  else
    echo "${!i} not known. Exiting."
    exit 1
  fi
  i=$((i+1))
done
ind=`indent $nesting_level "#"`
''')
    fp.write('echo "$ind Postscript {} begins."\n'.format(postscriptname))
    fp.write('echo "$ind Completing the task-'+r'${TASK}'+' build of {}"\n'.format(currpsf))
    vmd_instructions(fp,psfgen,logname=r'psfgen${TASK}.log',msg='generates psf={} pdb={}'.format(currpsf,currpdb))
    fp.write("cat {} | sed \'1,/#### BEGIN PATCHES/d;/#### END PATCHES/,$d\' > patches.inp\n".format(psfgen))
    outname=r'postnamd${TASK}-1'
    currcfg=r'run${TASK}-1.namd'
    currlog=r'run${TASK}-1.log'
    namd_instructions(fp,currcfg,currpsf,currpdb,outname,currlog,npe=npe,
                      numminsteps=nummin,numsteps=numsteps,seed=random.randint(0,10000),
                      template='vac.namd',temperature=temperature,msg='first relaxation',
                      stdparamfiles=StdParamFiles,localparamfiles=LocalParamFiles)
    namdbin='{}.coor'.format(outname)
    currpdb='{}.pdb'.format(outname)
    vmd_instructions(fp,r'$PSFGEN_BASEDIR/scripts/namdbin2pdb.tcl',logname=r'namdbin2pdb${TASK}-1.log',
                        args='{} {} {}'.format(currpsf,namdbin,'tmp.pdb'),msg='converting namdbin to pdb')
    fp.write('cat charmm_header.pdb tmp.pdb > {}\n'.format(currpdb))
    logname=r'ringp${TASK}-1.log'
    vmd_instructions(fp,r'$PSFGEN_BASEDIR/scripts/ringp.tcl',logname=logname,args='{} {}'.format(currpsf,currpdb),msg='checking for pierced rings')
    fp.write('npiercings=`grep -c pierces {}`\n'.format(logname))
    fp.write(r'if [[ $npiercings -gt 0 ]]; then'+'\n')
    fp.write(r'  echo "Error: There are $npiercings piercings in '+'{}"\n'.format(currpdb))
    fp.write('  grep pierces {}\n'.format(logname))
    fp.write('  echo "$ind Change your relaxation parameters and try again."\n')
    fp.write('  exit 1\n')
    fp.write('fi\n')
    if 'do_preclose_min_smd' in PostMod and PostMod['do_preclose_min_smd']:
        temperature_close=400
        target_numsteps=20000
        if 'preclose_params' in PostMod:
            p=PostMod['preclose_params']
            temperature_close=DefOrDict(p,'temperature_close',temperature_close)
            target_numsteps=DefOrDict(p,'target_numsteps',target_numsteps)
        fp.write('cat > close_these.inp << EOF\n')
        for l in sorted(Loops, key=lambda x: len(x.residues)):
            if (l.term and len(l.residues)>2):
                fp.write('{} {} {}\n'.format(l.replica_chainID,l.residues[-1].resseqnum,l.nextfragntermresid))
        fp.write('EOF\n')
        # measures to find the initial distances; generated fixed.pdb to fix the N atoms 
        logname=r'close${TASK}.log'
        args='{} {} close_these.inp fixed.pdb'.format(currpsf,currpdb)
        vmd_instructions(fp,r'$PSFGEN_BASEDIR/scripts/measure_bonds.tcl',logname=logname,args=args,msg='creating closing input')
        fp.write('if [ -f cv.inp ]; then rm cv.inp; fi\n')
        fp.write('touch cv.inp\n')
        fp.write('while IFS=" " read -r C L R B; do\n')
        fp.write(r'  cat $PSFGEN_BASEDIR/templates/cv-template.in | sed s/%C%/$C/g |')
        fp.write(r'  sed s/%NAME%/${C}${L}/g | sed s/%I%/$L/g | sed s/%J%/$R/g | sed s/%R0%/$B/g |')
        fp.write('  sed s/%TARGETNUMSTEPS%/{}/ >> cv.inp ;\n'.format(target_numsteps))
        fp.write('done < close_these.inp\n')
        outname=r'postnamd${TASK}-2'
        currcfg=r'run${TASK}-2.namd'
        currlog=r'run${TASK}-2.log'
        extras=['fixedatoms on','fixedatomsfile fixed.pdb','fixedatomscol B','colvars on','colvarsconfig cv.inp']
        namd_instructions(fp,currcfg,currpsf,currpdb,outname,currlog,npe=npe,
                      numminsteps=0,numsteps=int(1.5*target_numsteps),seed=random.randint(0,10000),
                      template='vac.namd',temperature=temperature_close,extras=extras,msg='closing',
                      stdparamfiles=StdParamFiles,localparamfiles=LocalParamFiles)
        namdbin='{}.coor'.format(outname)
        currpdb='{}.pdb'.format(outname)
        vmd_instructions(fp,r'$PSFGEN_BASEDIR/scripts/namdbin2pdb.tcl',logname=r'namdbin2pdb${TASK}-1.log',
                        args='{} {} {}'.format(currpsf,namdbin,'tmp.pdb'),msg='converting namdbin to pdb')
        fp.write('cat charmm_header.pdb tmp.pdb > {}\n'.format(currpdb))
        fp.write('cat > closure_patches.inp << EOF\n')
        for l in sorted(Loops, key=lambda x: len(x.residues)):
            if (l.term and len(l.residues)>2):
                #fp.write('# will try to close bond between {} and {} on chain {}...\n'.format(l.residues[-1].resseqnum,l.nextfragntermresid,l.replica_chainID))
                fp.write('patch HEAL {c}:{ll} {c}:{l} {c}:{r} {c}:{rr}\n'.format(c=l.replica_chainID,
                            ll=l.residues[-2].resseqnum,l=l.residues[-1].resseqnum,r=l.nextfragntermresid,rr=(l.nextfragntermresid+1)))
        fp.write('EOF\n')
        tfp=open('topologies.inp','w')
        CommonPSFGENheader(tfp,CTopo,LocTopo)
        tfp.close()
        fp.write('cat $PSFGEN_BASEDIR/scripts/loop_closure.tcl | sed "/#### LIGATION LIST STARTS/r the_closing_patches.inp"')
        fp.write(' | sed "/#### TOPOLOGY FILE LIST STARTS/r topologies.inp" > do_the_closures.tcl\n')
        newpsf='ligated.psf'
        newpdb='ligated.pdb'
        vmd_instructions(fp,'do_the_closures.tcl',args='{} {} {} {}'.format(currpsf,currpdb,newpsf,newpdb),logname=r'ligations${TASK}.log')
        currpsf=newpsf
        currpdb=newpdb
        currcfg=r'run${TASK}-3.namd'
        currlog=r'run${TASK}-3.log'
        outname=r'postnamd${TASK}-3'
        namd_instructions(fp,currcfg,currpsf,currpdb,outname,currlog,npe=npe,
                      numminsteps=nummin,numsteps=numsteps,seed=random.randint(0,10000),
                      template='vac.namd',temperature=temperature,msg='minimization of ligated peptide bonds',
                      stdparamfiles=StdParamFiles,localparamfiles=LocalParamFiles)
        namdbin='{}.coor'.format(outname)
        currpdb='{}.pdb'.format(outname)
        vmd_instructions(fp,r'$PSFGEN_BASEDIR/scripts/namdbin2pdb.tcl',logname=r'namdbin2pdb${TASK}-1.log',
                        args='{} {} {}'.format(currpsf,namdbin,'tmp.pdb'),msg='converting namdbin to pdb')
        fp.write('cat charmm_header.pdb tmp.pdb > {}\n'.format(currpdb))
 
        logname=r'ringp${TASK}-2.log'
        vmd_instructions(fp,r'$PSFGEN_BASEDIR/scripts/ringp.tcl',logname=logname,args='{} {}'.format(currpsf,currpdb),msg='checking for pierced rings')
        fp.write('npiercings=`grep -c pierces {}`\n'.format(logname))
        fp.write(r'if [[ $npiercings -gt 0 ]]; then'+'\n')
        fp.write(r'  echo "Error: There are $npiercings piercings in '+'{}"\n'.format(currpdb))
        fp.write('  grep pierces {}\n'.format(logname))
        fp.write('  echo "Change your relaxation parameters and try again."\n')
        fp.write('  exit 1\n')
        fp.write('fi\n')

    fp.write('echo {} {} {} > .tmpvar\n'.format(currpsf,currpdb,currcfg))
    fp.write('echo "$ind Postscript {} finishes."\n'.format(postscriptname))
    fp.write('exit 0')
    fp.close()
    os.system('chmod 744 {}'.format(postscriptname))

    
