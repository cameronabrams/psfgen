import sys
import operator
from datetime import date 
from atom import Atom,_PDBAtomNameDict_
from residue import Residue, _PDBResName123_, _pdb_glycans_, _pdb_ions_,_PDBResNameDict_
from chain import Chain
from segment import Segment, _seg_class_
from link import Link,LinkSet
''' 
    Parses experimental PDB to build input file for VMD/psfgen
    Cameron F Abrams
    cfa22@drexel.edu
'''
class Missing:
    def __init__(self,record_name,code,model,resname,chainID,resseqnum):
        self.record_name=record_name
        self.code=code
        self.model=model
        self.resname=resname
        self.chainID=chainID
        self.resseqnum=resseqnum
    def __str__(self):
        retstr='MISSING\n'+\
               '   model     {:s}\n'+\
               '   resname   {:s}\n'+\
               '   chainID   {:s}\n'+\
               '   resseqnum {:d}\n'
        return retstr.format(self.model,self.resname,self.chainID,self.resseqnum)
    def psfgen_residueline(self):
        return '     residue {} {} {}'.format(self.resname,self.resseqnum,self.chainID)

def read_missing(line):
    record_name=line[0:6]
    code=int(line[7:10])
    model=line[13:14]
    resname=line[15:18]
    chainID=line[19:20]
    resseqnum=int(line[21:26])
#    if (resseqnum<0):
#       fp.write('# negative resid for missing residue: {} {} {}'.format(resname,chainID,resseqnum))
    return(Missing(record_name,code,model,resname,chainID,resseqnum))    
def make_chains(Residues):
    C=[]
    for r in Residues:
        if C==[]:
            C.append(Chain(r))
        else:
            found=False
            for c in C:
                if c.chainID==r.chainID:
                    found=True
                    c.add_residue(r)
                    break
            if not found:
               C.append(Chain(r))
    for c in C:
        c.sort_residues()
    return C
def make_residues(Atoms,Missing):
    R=[]
    r=0
    for a in Atoms:
        if r==0:
            R.append(Residue(a=a))
            r=R[-1]
        else:
            if r.resseqnum==a.resseqnum and r.name == a.resname and r.chainID==a.chainID:
                r.add_atom(a=a)
            else:
                R.append(Residue(a=a))
                r=R[-1]
    for m in Missing:
        R.append(Residue(m=m))
    return R

class Mutation:
    def __init__(self,chainID,lr,ri,rr,label):
        self.chainID=chainID
        self.orig=lr
        self.resseqnum=ri
        self.new=rr
        self.label=label
    def __str__(self):
        return '{}-{}{}{}'.format(self.chainID,self.orig,self.resseqnum,self.new)
    def psfgen_segment_str(self):
        #return '    if {{ ${} == 1 }} {{\n        mutate {} {}\n    }}\n'.format(self.label,self.resseqnum,self.new)
        return '   mutate {} {}\n'.format(self.resseqnum,self.new)

def read_mutation_user(label):
    if not label[0].isalpha() or label[1]!='_':
        print('Poorly formed mutation spec: {}'.format(label))
        return 0
    chainID=label[0]
    lr=_PDBResName123_[label[2]]
    rr=_PDBResName123_[label[-1]]
    ri=int(label[3:-1])
    return Mutation(chainID,lr,ri,rr,label)

class Cleavage:
    def __init__(self,parent_chainID,parent_Cterm_resseqnum,daughter_chainID):
        self.parent_chainID=parent_chainID
        self.parent_Cterm_resseqnum=parent_Cterm_resseqnum
        self.daughter_chainID=daughter_chainID
    def __str__(self):
        return '{}{}-x-{}'.format(self.parent_chainID,self.parent_Cterm_resseqnum,self.daughter_chainID)

def read_cleavage_user(label):
    if not label[0].isalpha() or not label[-1].isalpha():
        print('Poorly formed cleavage spec: {}'.format(label))
        return 0
    parent_chainID=label[0]
    parent_Cterm_resseqnum=int(label[1:-1])
    daughter_chainID=label[-1]
    return Cleavage(parent_chainID,parent_Cterm_resseqnum,daughter_chainID)

class Molecule:
    def __init__(self,index,pdb):
        self.index=index
        self.pdb=pdb
        self.Atoms=[]
        self.Links=LinkSet()
        self.SSBonds=[]
        self.MissingRes=[]
        with open(pdb) as pdbfile:
            for line in pdbfile:
                if line[:4] == 'ATOM' or line[:6] == "HETATM":
                    self.Atoms.append(Atom(line))
                elif line[:4] == 'LINK':
                    self.Links.add_link(Link(line))
                elif line[:6] == 'SSBOND':
                    ss=read_ssbond(line)
                    self.SSBonds.append(ss)
                elif line[:6] == 'REMARK':
                    code=int(line[7:10])
                    test_int=line[20:26].strip()
                    if code==465 and (test_int.isdigit() or (len(test_int)>0 and  test_int[0]=='-')):
                        mr=read_missing(line)
                        self.MissingRes.append(mr)
        self.Residues=make_residues(self.Atoms,self.MissingRes)
        self.Chains=make_chains(self.Residues)
        
    def Cleave(self,Cleavages):
        for clv in Cleavages:
            clv_c=-1
            daughter_chainID_ok=True
            for c in self.Chains:
                if c.chainID==clv.parent_chainID:
                    clv_c=c
                if c.chainID==clv.daughter_chainID:
                    daugher_chainID_ok=False
            if clv_c!=-1 and daughter_chainID_ok:
                print('### before cleave:',clv_c)
                daughter=clv_c.Cleave(clv)
                self.Chains.append(daughter)
                for s in self.SSBonds:
                    if s.chainID1==clv_c.chainID and s.resseqnum1>clv.parent_Cterm_resseqnum:
                        s.chainID1=daughter.chainID
                    if s.chainID2==clv_c.chainID and s.resseqnum2>clv.parent_Cterm_resseqnum:
                        s.chainID2=daughter.chainID
                for l in self.Links:
                    if l.chainID1==clv_c.chainID and l.resseq1>clv.parent_Cterm_resseqnum:
                        l.chainID1=daughter.chainID
                    if l.chainID2==clv_c.chainID and l.resseq2>clv.parent_Cterm_resseqnum:
                        l.chainID2=daughter.chainID
                # to do -- links!   
                print('### after cleave:',clv_c,self.Chains[-1])
            else:
                print('### unable to cleave chain {} at position {} to generate {} {}'.format(clv_c.chainID,clv.parent_Cterm_resseqnum,clv.parent_chainID,clv.daughter_chainID))
    def __str__(self):
        return 'Molecule {} {}: {} chains, {} residues, {} atoms, {} links, {} ssbonds'.format(self.index,self.pdb,len(self.Chains),len(self.Residues),len(self.Atoms),len(self.Links.L),len(self.SSBonds)) 
    def residue_shift(self,chainID,resseqnumshift):
        found=False
        for c in self.Chains:
            if c.chainID==chainID:
                 found=True
                 break;
        if not found:
            print('### Could not apply shift to chain {}: no such chain in Molecule {}'.format(chainID,self.index))
            return -1
        for r in c.residues:
            r.residue_shift(resseqnumshift)
        for l in self.Links.L:
            if l.chainID1==chainID:
                l.resseqnum1+=resseqnumshift
            if l.chainID2==chainID:
                l.resseqnum2+=resseqnumshift            
        for ss in self.SSBonds:
            if ss.chainID1==chainID:
                ss.resseqnum1+=resseqnumshift
            if ss.chainID2==chainID:
                ss.resseqnum2+=resseqnumshift
        return 0

    def writepsfgeninput(self,fp,Mutations=[],topologies=[]):
        fp.write('if {![info exists PSFGEN_BASEDIR]} {\n'+\
              '    if {[info exists env(PSFGEN_BASEDIR]} {\n'+\
              '        set PSFGEN_BASEDIR $env(PSFGEN_BASEDIR)\n'+\
              '    } else {\n'+\
              '        set PSFGEN_BASEDIR $env(HOME)/research/psfgen\n'+\
              '    }\n'+\
              '}\n')
        fp.write('if {![info exists CHARMM_TOPPARDIR]} {\n'+\
              '    if {[info exists env(CHARMM_TOPPARDIR]} {\n'+\
              '        set TOPPARDIR $env(CHARMM_TOPPARDIR)\n'+\
              '    } else {\n'+\
              '        set TOPPARDIR $env(HOME)/charmm/toppar\n'+\
              '    }\n'+\
              '}\n')
        fp.write('source ${PSFGEN_BASEDIR}/src/loopmc.tcl\n')
        fp.write('source ${PSFGEN_BASEDIR}/scripts/vmdrc.tcl\n')
        fp.write('package require psfgen\n')
        for t in topologies:
            fp.write('topology $TOPPARDIR/{}\n'.format(t))
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

        fp.write('mol new {}\n'.format(self.pdb))

        for k,v in _PDBResNameDict_.items():
            fp.write('set RESDICT({}) {}\n'.format(k,v))
        for k,v in _PDBAtomNameDict_.items():
            fp.write('set ANAMEDICT({}) {}\n'.format(k,v))

        fp.write('set logid -1\n')

        Loops=[]
        for c in self.Chains:
            for s in c.Segments(Mutations=Mutations):
                stan,supp,coor,caco,loops=s.psfgen_segmentstanza()
                fp.write('### begin stanza for segment {}\n'.format(s.segname))
                fp.write(supp)
                fp.write(stan+'\n')
                fp.write(coor)
                fp.write(caco)
                if len(loops)>0:
                    Loops.extend(loops)
                fp.write('### end stanza for segment {}\n'.format(s.segname))

        for ss in self.SSBonds:
            fp.write(ss.psfgen_patchline())

        for l in self.Links.L:
            fp.write(l.psfgen_patchline())

        fp.write('guesscoord\n')
        fp.write('regenerate angles dihedrals\n')

        prefix=self.pdb[:]
        prefix=prefix.replace('.pdb','')
        fp.write('writepsf my_{}.psf\n'.format(prefix))
        fp.write('writepdb my_{}_raw.pdb\n'.format(prefix))
        return Loops

def WritePostMods(fp,pdb,center_protein,reorient_protein,reorselstr,do_loop_mc,Loops):
    prefix=pdb[:]
    prefix=prefix.replace('.pdb','')
    fp.write('mol delete top\n')
    fp.write('mol new my_{}.psf\n'.format(prefix))
    fp.write('set molid [molinfo top get id]\n')
    fp.write('mol addfile my_{}_raw.pdb\n'.format(prefix))
    if center_protein:
        fp.write('set a [atomselect top "all"]\n')
        fp.write('set or [measure center $a weight mass]\n')
        fp.write('$a moveby [vecscale -1 $or]\n')
        if reorient_protein:
            fp.write('set ca [measure center [atomselect top "protein and {}"] weight mass]\n'.format(reorselstr[0]))
            fp.write('set cb [measure center [atomselect top "protein and {}"] weight mass]\n'.format(reorselstr[1]))
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
    if do_loop_mc:
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

    fp.write('$a writepdb my_{}.pdb\n'.format(prefix))

if __name__=='__main__':

    print('### cfapdbparser {}'.format(date.today()))
    i=1
    Molecules=[]
    Mut=[]
    Clv=[]
    do_loop_mc=False
    center_protein=True
    reorient_protein=False
    reorselstr=[]
    psfgen='mkpsf.tcl'
    topo=['top_all36_prot.rtf','top_all36_carb_namd_cfa.rtf','stream/carb/toppar_all36_carb_glycopeptide.str','toppar_water_ions_namd_nonbfixes.str']
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
                mut=read_mutation_user(sys.argv[i])
                if mut!=0:
                    Mut.append(mut)
                i+=1
            if i<len(sys.argv) and sys.argv[i][0]=='-':
                i-=1
        elif sys.argv[i]=='-cleave':
            i+=1
            while i<len(sys.argv) and sys.argv[i][0]!='-':
                clv=read_cleavage_user(sys.argv[i])
                if clv!=0:
                    Clv.append(clv)
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
            do_loop_mc=True
        elif sys.argv[i]=='-no_center':
            center_protein=False
        elif sys.argv[i]=='-reorient_protein':
            reorient_protein=True
            i+=1
            while i<len(sys.argv) and sys.argv[i][0]!='-':
                reorselstr.append(sys.argv[i])
                i+=1
            if sys.argv[i][0]=='-':
                i-=1
        elif sys.argv[i]=='-psfgen':
            i+=1
            psfgen=sys.argv[i]
        i+=1

    print('### psfgen input to be created {}'.format(psfgen))

    if reorient_protein:
        center_protein=True
        if len(reorient_selstr)<2:
            print('Error: must specify two atomselections to define local-z axis for reorientation')
            print('       disabling reorientation.')
            reorient_protein=False

    Base=Molecules[0]
    
    # do stuff to Base molecule using stuff from others
    if len(Clv)>0:
        Base.Cleave(Clv)
        print('### after cleavages:')
        for c in Base.Chains:
            print(c)

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
    Loops=Base.writepsfgeninput(psfgen_fp,Mut,topo)
    WritePostMods(psfgen_fp,Base.pdb,center_protein,reorient_protein,reorselstr,do_loop_mc,Loops)
    psfgen_fp.write('exit\n')
    psfgen_fp.write('### thank you for using cfapdbparser!\n')
    print('### next: vmd -dispdev text -e {}'.format(psfgen))
    
