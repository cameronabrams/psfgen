import sys
import operator
from datetime import date 

''' Parses experimental PDB to 
       1. Extract loops of missing structure from REMARK 465 lines in a PDB
       2. Extract N-linked NAGs (more to come)
       3. Extract disulfides

    stdout can be redirected to a file that can be sourced by a mkpsf script
    after PDB is read in but before the 'guesscoords' statement

    Cameron F Abrams
    cfa22@drexel.edu
'''
_pdb_ions_=['LIT','SOD','MG','POT','CAL','RUB','CES','BAR','ZN','CAD','CL']
_pdb_glycans_=['BMA','FUC','GAL','MAN','NAG','SIA']

_PDBResName123_={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY',
               'H':'HSE','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER',
               'T':'THR','W':'TRP','Y':'TYR','V':'VAL'}

_seg_class_={'HOH':'WATER'}
_seg_class_.update({k:'ION' for k in _pdb_ions_})
_seg_class_.update({k:'GLYCAN' for k in _pdb_glycans_})
_seg_class_.update({k:'PROTEIN' for k in _PDBResName123_.values()})
_seg_class_['HIS']='PROTEIN'
_segname_suffix_character_={'PROTEIN':'','ION':'I','WATER':'WX','GLYCAN':'S'}
_PDBResNameDict_={'HIS':'HSE','ZN':'ZN2','HOH':'TIP3','CL':'CLA','NAG':'BGNA','MAN':'AMAN','BMA':'BMAN','FUC':'AFUC','GAL':'BGAL','SIA':'ANE5AC'}
_PDBAtNameDict_={'CL':'CLA'}



class Atom:
   def __init__(self,record_name,serial,name,altloc,resname,chainID,resseqnum,insertion,x,y,z,occ,beta,elem,charge):
      self.record_name=record_name
      self.serial=serial
      self.name=name
      self.altloc=altloc
      self.resname=resname
      self.chainID=chainID
      self.resseqnum=resseqnum
      self.insertion=insertion
      self.x=x
      self.y=y
      self.z=z
      self.occ=occ
      self.beta=beta
      self.elem=elem
      self.charge=charge
   def __str__(self):
      retstr='{}\n'+\
             '  serial    {:d}\n'+\
             '  name      {:s}\n'+\
             '  altloc    {:s}\n'+\
             '  resname   {:s}\n'+\
             '  chainID   {:s}\n'+\
             '  resseqnum {:d}\n'+\
             '  insertion {:s}\n'+\
             '  x         {:.3f}\n'+\
             '  y         {:.3f}\n'+\
             '  z         {:.3f}\n'+\
             '  occ       {:.2f}\n'+\
             '  beta      {:.2f}\n'+\
             '  elem      {:s}\n'+\
             '  charge    {:s}\n'
      return retstr.format(self.record_name,self.serial,self.name,self.altloc,self.resname,self.chainID,self.resseqnum,self.insertion,self.x,self.y,self.z,self.occ,self.beta,self.elem,self.charge)

class Link:
    def __init__(self,record_name,name1,altloc1,resname1,chainID1,resseq1,icode1,name2,altloc2,resname2,chainID2,resseq2,icode2,sym1,sym2,link_distance):
       self.record_name=record_name.strip()
       self.name1=name1.strip()
       self.altloc1=altloc1.strip()
       self.resname1=resname1.strip()
       self.chainID1=chainID1.strip()
       self.resseq1=resseq1
       self.icode1=icode1.strip()
       self.name2=name2.strip()
       self.altloc2=altloc2.strip()
       self.resname2=resname2.strip()
       self.chainID2=chainID2.strip()
       self.resseq2=resseq2
       self.icode2=icode2.strip()
       self.sym1=sym1.strip()
       self.sym2=sym2.strip()
       self.link_distance=link_distance
    def __str__(self):
        retstr='{}\n'+\
                '   name1         {:s}\n'+\
                '   altloc1       {:s}\n'+\
                '   resname1      {:s}\n'+\
                '   chainID1      {:s}\n'+\
                '   resseq1       {:d}\n'+\
                '   icode1        {:s}\n'+\
                '   name2         {:s}\n'+\
                '   altloc2       {:s}\n'+\
                '   resname2      {:s}\n'+\
                '   chainID2      {:s}\n'+\
                '   resseq2       {:d}\n'+\
                '   icode2        {:s}\n'+\
                '   sym1          {:s}\n'+\
                '   sym2          {:s}\n'+\
                '   link_distance {:.3f}'
        return retstr.format(self.record_name,self.name1,self.altloc1,self.resname1,self.chainID1,self.resseq1,self.icode1,self.name2,self.altloc2,self.resname2,self.chainID2,self.resseq2,self.icode2,self.sym1,self.sym2,self.link_distance)
    def psfgen_patchline(self):
        if self.resname1=='ASN' and _seg_class_[self.resname2]=='GLYCAN':
            return 'patch NGLB {}:{} {}S:{}'.format(self.chainID1,self.resseq1,self.chainID2,self.resseq2)
        else:
            # for a glycan-glycan patch, the C1 atom is always on the i-residue
            if self.name2=='C1' and _seg_class_[self.resname1]=='GLYCAN':
                cmdi='[axeq {} 0 {} {} {}]'.format(self.resseq2,self.chainID2,self.name2,self.resseq1)
                cmdj='[axeq {} 0 {} {} {}]'.format(self.resseq1,self.chainID1,self.name1,-1)
                
                return 'patch 1{:1s}{}{} {}S:{} {}S:{}'.format(self.name1[1], cmdi,cmdj,self.chainID2,self.resseq2,self.chainID1,self.resseq1)
            elif self.name1=='C1' and _seg_class_[self.resname2]=='GLYCAN':
                return 'patch 1{:1s}xx {}S:{} {}S:{}'.format(self.name2[1], self.chainID1,self.resseq1,self.chainID2,self.resseq2)
            else:
                return 'patch unknown for '+str(self)

class SSBond:
    def __init__(self,record_name,resname1,chainID1,resseqnum1,icode1,resname2,chainID2,resseqnum2,icode2,sym1,sym2,length):
        self.record_name=record_name.strip()
        self.resname1=resname1.strip()
        self.chainID1=chainID1.strip()
        self.resseqnum1=resseqnum1
        self.icode1=icode1.strip()
        self.resname2=resname2.strip()
        self.chainID2=chainID2.strip()
        self.resseqnum2=resseqnum2
        self.icode2=icode2.strip()
        self.sym1=sym1.strip()
        self.sym2=sym2.strip()
        self.length=length
    def __str__(self):
        retstr='{}\n'+\
               '  resname1    {:s}\n'+\
               '  chainID1    {:s}\n'+\
               '  resseqnum1  {:d}\n'+\
               '  icode1      {:s}\n'+\
               '  resname2    {:s}\n'+\
               '  chainID2    {:s}\n'+\
               '  resseqnum2  {:d}\n'+\
               '  icode2      {:s}\n'+\
               '  sym1        {:s}\n'+\
               '  sym2        {:s}\n'+\
               '  length      {:.3f}\n'
        return retstr.format(self.record_name,self.resname1,self.chainID1,self.resseqnum1,self.icode1,self.resname2,self.chainID2,self.resseqnum2,self.icode2,self.sym1,self.sym2,self.length)
    def psfgen_patchline(self):
       return 'patch DISU {}:{} {}:{}'.format(self.chainID1,self.resseqnum1,self.chainID2,self.resseqnum2)

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

def psfgen_write_pdb(st,c,l,r,p):
    if st=='PROTEIN':
        retstr='[atomselect top "chain {} and resid {} to {}"] writepdb {}\n'.format(c,l,r,p)
    else:
        retstr='# Changing resnames and atom names in {}\n'.format(p)
        retstr+='set myseg [atomselect top "chain {} and resid {} to {}"]\n'.format(c,l,r,p)
        retstr+='set sav_nm [$myseg get resname]\n'
        retstr+='set new_nm [list]\n'
        retstr+=r'foreach r $sav_nm {'+'\n'
        retstr+='   lappend new_nm $RESDICT($r)\n'
        retstr+='}\n'
        retstr+='$myseg set resname $new_nm\n'
        retstr+='set new_nm [list]\n'
        retstr+='set sav_nm [$myseg get name]\n'
        retstr+=r'foreach r $sav_nm {'+'\n'
        retstr+=r'   if { [ info exists ANAMEDICT($r) ] } {'+'\n'
        retstr+='      lappend new_nm $ANAMEDICT($r)\n'
        retstr+=r'   } else {'+'\n'
        retstr+='      lappend new_nm $r\n'
        retstr+='   }\n'
        retstr+='}\n'
        retstr+='$myseg set name $new_nm\n'
        if st=='WATER':
            retstr+='$myseg set name OH2\n'
        retstr+='$myseg writepdb {}\n'.format(p)

    return retstr

class Run:
    def __init__(self,chainID,resseqnum1,resseqnum2):
        self.chainID=chainID
        self.resseqnum1=resseqnum1
        self.resseqnum2=resseqnum2
    def pdb_str(self):
        return '{}_{}_to_{}.pdb'.format(self.chainID,self.resseqnum1,self.resseqnum2)
    def __str__(self):
        return 'RUN: {} {} to {}'.format(self.chainID,self.resseqnum1,self.resseqnum2)
class Loop:
    def __init__(self,chainID,resseqnum0,r):
        self.chainID=chainID
        self.resseqnum0=resseqnum0
        self.residues=[r]
        self.terminated=False
    def add_residue(self,r):
        self.residues.append(r)
    def __str__(self):
        return 'LOOP chainID {} r0 {} r1 {} to r2 {}'.format(self.chainID,self.resseqnum0,self.residues[0].resseqnum,self.residues[-1].resseqnum)
    def caco_str(self):
        return 'coord {} {} N [cacoIn_nOut {} {} 0]\n'.format(self.chainID,self.residues[0].resseqnum,self.resseqnum0,self.chainID)

class Segment:
    def __init__(self,r):
        self.segname=r.chainID+_segname_suffix_character_[_seg_class_[r.name]]
        self.segtype=_seg_class_[r.name]
        self.residues=[r]
        self.mutations=[]
    def add_residue(self,r):
        self.residues.append(r)
    def psfgen_segmentstanza(self):  # good job, now this only works for proteins!!
        pdbs=[]
        Loops=[]
        Runs=[]
        cacostr=''
        suppstr=''
        retstr='segment {} {{\n'.format(self.segname)
        inloop=False
        last_resseqnum=-1
#        for r in self.residues:
#            print(r)
        for r in self.residues:
            if len(r.atoms)>0:
                if last_resseqnum==-1:
#                    print('appending',r)
                    Runs.append(Run(r.chainID,r.resseqnum,-1))
                elif inloop:
                    if len(Loops)>0:
                        Loops[-1].terminated=True
                    inloop=False
                    Runs.append(Run(r.chainID,r.resseqnum,-1))
                else:
#                    print('2passing',r)
                    pass
            elif len(r.atoms)==0:
                if last_resseqnum==-1:
                    pass # Loops.append(Loop(r.chainID,-1,r))                    
                elif not inloop:
                    Runs[-1].resseqnum2=last_resseqnum
                    Loops.append(Loop(r.chainID,last_resseqnum,r))
                else:
#                    print(last_resseqnum)
                    if len(Loops)>0:
                        Loops[-1].add_residue(r)
                    else:
                        pass
                inloop=True
            else:
#                print('passing', r)
                pass
            last_resseqnum=r.resseqnum
#        print(self.segname,len(Runs))
        if len(Runs)==1:
            Runs[0].resseqnum2=last_resseqnum
#            print('END OF LOOP: segname {} r.resseqnum {} {} last_resseqnum {}'.format(self.segname,r.resseqnum,r.name,last_resseqnum))
        if len(Loops)>0:
            for r,l in zip(Runs,Loops):
                pdbs.append(r.pdb_str())
                suppstr+=psfgen_write_pdb(self.segtype,r.chainID,r.resseqnum1,r.resseqnum2,pdbs[-1])
                retstr+='   pdb {}\n'.format(pdbs[-1])
                if l.terminated==True:
                    for rr in l.residues:
                        if rr.name in _PDBResNameDict_:
                            nm=_PDBResNameDict_[rr.name]
                        else:
                            nm=rr.name
                        retstr+='   residue {} {} {}\n'.format(rr.resseqnum,nm,rr.chainID)
            for m in self.mutations:
                retstr+=m.psfgen_segment_str()
        else:
            r=Runs[0]
            pdbs.append(r.pdb_str())
            suppstr+=psfgen_write_pdb(self.segtype,r.chainID,r.resseqnum1,r.resseqnum2,pdbs[-1])
            retstr+='   pdb {}\n'.format(pdbs[-1])

        retstr+='}'
        coordstr=''
        for p in pdbs:
            coordstr+='coordpdb {} {}\n'.format(p,s.segname)
        if len(Loops)>0:
            r0=self.residues[0]
            for l in Loops:
                if l.terminated:
                    cacostr+=l.caco_str()
        return retstr,suppstr,coordstr,cacostr,Loops

def read_atom(line):
# 1 -  6        Record name   "ATOM  "
    record_name=line[0:6]
# 7 - 11        Integer       serial       Atom  serial number.
    serial=int(line[6:11])
#13 - 16        Atom          name         Atom name.
    name=line[12:16]
#17             Character     altLoc       Alternate location indicator.
    altloc=line[16:17]
#18 - 20        Residue name  resName      Residue name.
    resname=line[17:20]
#22             Character     chainID      Chain identifier.
    chainID=line[21:22]
#23 - 26        Integer       resSeq       Residue sequence number.
    resseqnum=int(line[22:26])
#27             AChar         iCode        Code for insertion of residues.
    insertion=line[26:27]
#31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    x=float(line[30:38])
#39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    y=float(line[38:46])
#47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    z=float(line[46:54])
#55 - 60        Real(6.2)     occupancy    Occupancy.
    occ=float(line[54:60])
#61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    beta=float(line[60:66])
#77 - 78        LString(2)    element      Element symbol, right-justified.
    elem=line[76:78]
#79 - 80        LString(2)    charge       Charge  on the atom.
    charge=line[78:80]
    new_at=Atom(record_name,serial,name,altloc,resname,chainID,resseqnum,insertion,x,y,z,occ,beta,elem,charge)
    return new_at

def read_link(line):
# 1 -  6         Record name    "LINK  "
    record_name=line[0:6]
#13 - 16         Atom           name1           Atom name.
    name1=line[12:16]
#17              Character      altLoc1         Alternate location indicator.
    altloc1=line[16:17]
#18 - 20         Residue name   resName1        Residue  name.
    resname1=line[17:20]
#22              Character      chainID1        Chain identifier.
    chainID1=line[21:22]
#23 - 26         Integer        resSeq1         Residue sequence number.
    resseq1=int(line[22:26])
#27              AChar          iCode1          Insertion code.
    icode1=line[26:27]
#43 - 46         Atom           name2           Atom name.
    name2=line[42:46]
#47              Character      altLoc2         Alternate location indicator.
    altloc2=line[46:47]
#48 - 50         Residue name   resName2        Residue name.
    resname2=line[47:50]
#52              Character      chainID2        Chain identifier.
    chainID2=line[51:52]
#53 - 56         Integer        resSeq2         Residue sequence number.
    resseq2=int(line[52:56])
#57              AChar          iCode2          Insertion code.
    icode2=line[56:57]  
#60 - 65         SymOP          sym1            Symmetry operator atom 1.
    sym1=line[59:65]
#67 - 72         SymOP          sym2            Symmetry operator atom 2.
    sym2=line[66:72]
#74 – 78         Real(5.2)      Length          Link distance
    link_distance=float(line[73:78])
    return(Link(record_name,name1,altloc1,resname1,chainID1,resseq1,icode1,name2,altloc2,resname2,chainID2,resseq2,icode2,sym1,sym2,link_distance))

def read_ssbond(line):
# 1 -  6        Record name    "SSBOND"
    record_name=line[0:6]
# 8 - 10        Integer        serNum           Serial number.
    serial_number=int(line[7:10])
# 12 - 14        LString(3)     "CYS"            Residue name.
    resname1=line[11:14]
# 16             Character      chainID1         Chain identifier.
    chainID1=line[15:16]
# 18 - 21        Integer        seqNum1          Residue sequence number.
    resseqnum1=int(line[17:21])
# 22             AChar          icode1           Insertion code.
    icode1=line[21:22]
# 26 - 28        LString(3)     "CYS"            Residue name.
    resname2=line[25:28]
# 30             Character      chainID2         Chain identifier.
    chainID2=line[29:30]
# 32 - 35        Integer        seqNum2          Residue sequence number.
    resseqnum2=int(line[31:35])
# 36             AChar          icode2           Insertion code.
    icode2=line[35:36]
# 60 - 65        SymOP          sym1             Symmetry operator for residue 1.
    sym1=line[59:65]
# 67 - 72        SymOP          sym2             Symmetry operator for residue 2.
    sym2=line[66:72]
# 74 – 78        Real(5.2)      Length           Disulfide bond distance
    length=float(line[73:78])

    return(SSBond(record_name,resname1,chainID1,resseqnum1,icode1,resname2,chainID2,resseqnum2,icode2,sym1,sym2,length))

def read_missing(line):
    record_name=line[0:6]
    code=int(line[7:10])
    model=line[13:14]
    resname=line[15:18]
    chainID=line[19:20]
    resseqnum=int(line[22:26])
    return(Missing(record_name,code,model,resname,chainID,resseqnum))    
class Chain:
    def __init__(self,r):
        self.chainID=r.chainID
        self.residues=[r]
    def add_residue(self,r):
        if self.chainID==r.chainID:
            self.residues.append(r)
    def sort_residues(self):
        sorted_residues=sorted(self.residues,key=operator.attrgetter('resseqnum'))
        self.residues=sorted_residues
    def __str__(self):
        resid=[]
        for r in self.residues:
            resid.append(r.resseqnum)
        return 'CHAIN {} {:d} - {:d}'.format(self.chainID,min(resid),max(resid))
    def Segments(self,Mutations=0):
        S=[]
        for r in self.residues:
            if S==[]:
                S.append(Segment(r))
            else:
                for s in S:
                    #print('looking {} {}'.format(s.segtype,_seg_class_[r.name]))
                    found=False
                    if s.segtype==_seg_class_[r.name]:
                        #print('adding')
                        s.add_residue(r)
                        found=True
                        break
                if not found:
                    #print('creating')
                    S.append(Segment(r))
        if Mutations != 0:
            for m in Mutations:
                found=False
                if m.chainID==self.chainID:
                    for s in S:
                        if s.segtype=='PROTEIN':
                            for r in s.residues:
                                if r.resseqnum==m.resseqnum:
                                    found=True
                                    s.mutations.append(m)
                                    break
                                else:
                                    pass
                        else:
                            pass
                else:
                    pass
        else:
            pass

        return S

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

class Residue:
    def __init__(self,a=-1,m=0):
        if a!=-1:
            self.resseqnum=a.resseqnum
            self.name=a.resname
            self.chainID=a.chainID
            self.atoms=[a]
        else:
            if m!=0:
                self.resseqnum=m.resseqnum
                self.name=m.resname
                self.chainID=m.chainID
                self.atoms=[]
            else:
                print('ERROR: bad residue construction')

    def add_atom(self,a):
        if self.resseqnum==a.resseqnum and self.name==a.resname and self.chainID==a.chainID:
            self.atoms.append(a)
    def __str__(self):
        if len(self.atoms)==0:
            atstr='MISSING'
        else:
            atser=[]
            for a in self.atoms:
                atser.append(a.serial)
            atstr='{:d} - {:d}'.format(min(atser),max(atser))
        return 'RESIDUE {} {} {:d} {}'.format(self.chainID,self.name,self.resseqnum,atstr)

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

if __name__=='__main__':

    i=1
    pdb='null.pdb'
    Mutations=[]
    do_loop_mc=False
    center_protein=True
    reorient_protein=False
    reorselstr=[]
    topologies=['top_all36_prot.rtf','top_all36_carb_namd_cfa.rtf','stream/carb/toppar_all36_carb_glycopeptide.str','toppar_water_ions_namd_nonbfixes.str']
    while i<len(sys.argv):
        if sys.argv[i]=='-pdb':
            i+=1
            pdb=sys.argv[i]
        elif sys.argv[i]=='-mut':
            i+=1
            while i<len(sys.argv) and sys.argv[i][0]!='-':
                mut=read_mutation_user(sys.argv[i])
                if mut!=0:
                    Mutations.append(mut)
                i+=1
        elif sys.argv[i]=='-top':
            i+=1
            while i<len(sys.argv) and sys.argv[i][0]!='-':
                topologies.append(sys.argv[i])
                i+=1
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
        i+=1

    if reorient_protein:
        center_protein=True
        if len(reorient_selstr)<2:
            print('Error: must specify two atomselections to define local-z axis for reorientation')
            print('       disabling reorientation.')
            reorient_protein=False

    Atoms=[]
    Links=[]
    SSBonds=[]
    MissingRes=[]
    with open(pdb) as pdbfile:
        for line in pdbfile:
            if line[:4] == 'ATOM' or line[:6] == "HETATM":
                at=read_atom(line)
#            print(at)
                Atoms.append(at)
            elif line[:4] == 'LINK':
                ln=read_link(line)
#            print(ln.psfgen_patchline())
                Links.append(ln)
            elif line[:6] == 'SSBOND':
                ss=read_ssbond(line)
#            print(ss.psfgen_patchline())
                SSBonds.append(ss)
            elif line[:6] == 'REMARK':
                code=int(line[7:10])
                test_int=line[20:26].strip()
            #print(code,test_int)
                if code==465 and test_int.isdigit():
                    mr=read_missing(line)
#                print(mr.psfgen_residueline())
                    MissingRes.append(mr)

    Residues=make_residues(Atoms,MissingRes)
#    for r in Residues:
#       print(r)
    Chains=make_chains(Residues)

    # psfgen input begins here
    print('### cfapdbparser {}'.format(pdb))
    print('### {}'.format(date.today()))
    print('if {![info exists PSFGEN_BASEDIR]} {\n'+\
          '    if {[info exists env(PSFGEN_BASEDIR]} {\n'+\
          '        set PSFGEN_BASEDIR $env(PSFGEN_BASEDIR)\n'+\
          '    } else {\n'+\
          '        set PSFGEN_BASEDIR $env(HOME)/research/psfgen\n'+\
          '    }\n'+\
          '}\n')
    print('if {![info exists CHARMM_TOPPARDIR]} {\n'+\
          '    if {[info exists env(CHARMM_TOPPARDIR]} {\n'+\
          '        set TOPPARDIR $env(CHARMM_TOPPARDIR)\n'+\
          '    } else {\n'+\
          '        set TOPPARDIR $env(HOME)/charmm/toppar\n'+\
          '    }\n'+\
          '}\n')
    print(r'source ${PSFGEN_BASEDIR}/src/loopmc.tcl')
    print(r'source ${PSFGEN_BASEDIR}/scripts/vmdrc.tcl')
    print('package require psfgen')
    for t in topologies:
        print('topology $TOPPARDIR/{}'.format(t))


    print('pdbalias residue HIS HSD')
    print('pdbalias atom ILE CD1 CD')
    print('pdbalias residue NAG BGNA')
    print('pdbalias atom BGNA C7 C')
    print('pdbalias atom BGNA O7 O')
    print('pdbalias atom BGNA C8 CT')
    print('pdbalias atom BGNA N2 N')
    print('pdbalias residue SIA ANE5')
    print('pdbalias atom ANE5 C10 C')
    print('pdbalias atom ANE5 C11 CT')
    print('pdbalias atom ANE5 N5 N')
    print('pdbalias atom ANE5 O1A O11')
    print('pdbalias atom ANE5 O1B O12')
    print('pdbalias atom ANE5 O10 O')

    print('mol new {}'.format(pdb))

    for k,v in _PDBResNameDict_.items():
        print('set RESDICT({}) {}'.format(k,v))
    for k,v in _PDBAtNameDict_.items():
        print('set ANAMEDICT({}) {}'.format(k,v))

    print('set logid -1')

    Loops=[]
    for c in Chains:
        for s in c.Segments(Mutations=Mutations):
            stan,supp,coor,caco,loops=s.psfgen_segmentstanza()
            print('### begin stanza for segment {}'.format(s.segname))
            print(supp,end='')
            print(stan)
            print(coor,end='')
            print(caco,end='')
            if len(loops)>0:
                Loops.extend(loops)
            print('### end stanza for segment {}'.format(s.segname))

    for ss in SSBonds:
        print(ss.psfgen_patchline())

    for l in Links:
        print(l.psfgen_patchline())

    print('guesscoord')
    print('regenerate angles dihedrals')

    prefix=pdb[:]
    prefix=prefix.replace('.pdb','')
    print('writepsf my_{}.psf'.format(prefix))
    print('writepdb my_{}_raw.pdb'.format(prefix))
    print('mol delete top')
    print('mol new my_{}.psf'.format(prefix))
    print('set molid [molinfo top get id]')
    print('mol addfile my_{}_raw.pdb'.format(prefix))
    if center_protein:
        print('set a [atomselect top "all"]')
        print('set or [measure center $a weight mass]')
        print('$a moveby [vecscale -1 $or]')
        if reorient_protein:
            print('set ca [measure center [atomselect top "protein and {}"] weight mass]'.format(reorselstr[0]))
            print('set cb [measure center [atomselect top "protein and {}"] weight mass]'.format(reorselstr[1]))
            print('set pi 3.415928')
            print('set dv [vecsub $ca $cb]')
            print('set d [veclength $dv]')
            print('set cp [expr [lindex $dv 0]/$d]')
            print('set sp [expr [lindex $dv 1]/$d]')
            print('set p [expr acos($cp)]')
            print('if {[expr $sp < 0.0]} {')
            print('  set p [expr 2*$pi-$p]')
            print('}')
            print('set ct [expr [lindex $dv 2]/$d]')
            print('set t [expr acos($ct)]')
            print('$a move [transaxis z [expr -1 * $p] rad]')
            print('$a move [transaxis y [expr -1 * $t] rad]')
    if do_loop_mc:
        print('set loops {')
        for l in Loops:
            if l.terminated:
                print('{{ {} {} {} }}'.format(l.chainID,l.residues[0].resseqnum,l.residues[-1].resseqnum))
        print('           }')
        # create loops list { { }, { }, ...}
        print('set nc 1000')
        print('set rcut 3.0')
        print('set r0 1.5')
        print('set temperature 3.0')
        print('set k 10.0')
        print('set bg [atomselect $molid "noh"]')
        print('set loopindex 0')
        print('set nloops [llength $loops]')
        print('foreach l $loops {')
        print('   set chain [lindex $l 0]')
        print('   puts "Relaxing loop $loopindex out of $nloops"')
        print('   set residueList [[atomselect $molid "chain $chain and resid [lindex $l 1] to [lindex $l 2] and name CA"] get residue]')
        print('   do_loop_mc $residueList $chain $molid $k $r0 $bg $rcut $nc $temperature [irand_dom 1000 9999] $logid')
        print('   set loopindex [expr $loopindex + 1]')
        print('}')

    print('writepdb my_{}.pdb'.format(prefix))

    print('# {} links read.'.format(len(Links)))
    print('# {} ssbonds read.'.format(len(SSBonds)))
    print('# {} atoms read.'.format(len(Atoms)))
    print('# {} missing residues read.'.format(len(MissingRes)))
    print('# {} user-specfied mutations provided.'.format(len(Mutations)))
