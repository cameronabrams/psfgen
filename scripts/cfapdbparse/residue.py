_pdb_ions_=['LIT','SOD','MG','POT','CAL','RUB','CES','BAR','ZN','CAD','CL']
_pdb_glycans_=['BMA','FUC','GAL','MAN','NAG','SIA','ANE5']
_pdb_ligands_=['EIC','VCG']
_PDBResName123_={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY',
               'H':'HSE','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER',
               'T':'THR','W':'TRP','Y':'TYR','V':'VAL'}
_ResNameDict_PDB_to_CHARMM_={'HIS':'HSE','ZN':'ZN2','HOH':'TIP3','CL':'CLA',
'NAG':'BGNA','MAN':'AMAN','BMA':'BMAN','FUC':'AFUC','GAL':'BGAL','SIA':'ANE5AC','ANE5':'ANE5AC',
'EIC':'LIN','VCG':'VCG'}
#_ResNameDict_PDB_to_CHARMM_={'HIS':'HSE','ZN':'ZN2','HOH':'TIP3','CL':'CLA','NAG':'BGLC','MAN':'AMAN','BMA':'BMAN','FUC':'AFUC','GAL':'BGAL','SIA':'ANE5AC', 'ANE5':'ANE5AC'}
_ResNameDict_CHARMM_to_PDB_={v:k  for k,v in _ResNameDict_PDB_to_CHARMM_.items()}

def ResnameCharmify(nm):
    return _ResNameDict_PDB_to_CHARMM_[nm] if nm in _ResNameDict_PDB_to_CHARMM_ else nm

class Residue:
    def __init__(self,a=None,m=None,molid='*'):
        ''' initializing with an atom '''
        if a!=None:
            self.resseqnum=a.resseqnum
            self.insertion=a.insertion
            self.name=a.resname
            self.chainID=a.chainID
#            self.source_chainID=a.chainID
            self.atoms=[a]
            self.up=[]
            self.down=[]
            self.uplink=[]
            self.downlink=[]
        elif m!=None:
            ''' initializing as a missing residue '''
            self.resseqnum=m.resseqnum
            self.name=m.resname
            self.insertion=m.insertion
            self.chainID=m.chainID
#            self.source_chainID=m.chainID
            self.biomt='*'
            self.atoms=[]
            self.up=[]
            self.down=[]
            self.uplink=[]
            self.downlink=[]
        else:
            self.resseqnum=0
            self.insertion=' '
            self.name='UNK'
            self.chainID='UNK'
#            self.source_chainID
            self.atoms=[]
            self.up=[]
            self.down=[]
            self.uplink=[]
            self.downlink=[]
        self.resseqnumi=f'{self.resseqnum}{self.insertion}'
    def __lt__(self,other):
        if self.resseqnum<other.resseqnum:
            return True
        elif self.resseqnum==other.resseqnum:
            if self.insertion==None and other.insertion==None:
                return False
            elif (self.insertion=='' or self.insertion==' ' or self.insertion==None) and other.insertion.isalpha():
                return True
            elif self.insertion.isalpha() and other.insertion.isalpha():
                return ord(self.insertion)<ord(other.insertion)
            else:
                return False
    def add_atom(self,a):
        if self.resseqnum==a.resseqnum and self.name==a.resname and self.chainID==a.chainID:
            self.atoms.append(a)
    def residue_shift(self,resseqnumshift):
        self.resseqnum+=resseqnumshift
        for a in self.Atoms:
            a.resseqnum+=resseqnumshift
        pass
    def set_chainID(self,chainID):
        self.chainID=chainID
        for a in self.atoms:
            a.chainID=chainID
    def linkTo(self,other,link):
        self.down.append(other)
        self.downlink.append(link)
        other.up.append(self)
        other.uplink.append(link)
    def unlink(self,other,link):
        self.down.remove(other)
        self.downlink.remove(link)
        other.up.remove(self)
        other.uplink.remove(link)
    def ri(self):
        ins0='' if self.insertion==' ' else self.insertion
        return f'{self.resseqnum}{ins0}'
    def set_connections(self):
        pass
    def __str__(self):
        return '{}-{}{}'.format(self.chainID,self.name,self.resseqnum)
    def printshort(self):
        return str(self)
    def str_full(self):
        if len(self.atoms)==0:
            atstr='MISSING'
        else:
            atser=[]
            for a in self.atoms:
                atser.append(a.serial)
            atstr='{:d} - {:d}'.format(min(atser),max(atser))
        return '[{}] RESIDUE {} {} {:d} {}'.format(self.biomt,self.chainID,self.name,self.resseqnum,atstr)

    def get_down_group(self):
        res=[]
        for d in self.down:
            res.append(d)
            res.extend(d.get_down_group())
        return res

def get_residue(R,chainID,resseqnum,insertion=' '):
    candidates=[]
    #print(f'Searching list of {len(R)} residues for c {chainID} r {resseqnum}')
    for r in R:
        if r.chainID==chainID and r.resseqnum==resseqnum:
            candidates.append(r)
    if len(candidates)==1:
        return candidates[0]
    else:
        for r in candidates:
            if r.insertion==insertion:
                return r
    return None

def get_atom(R,chainID,resseqnum,atname,insertion=' '):
#    print('get_atom() searching for {} in resid {} chain {}'.format(atname,resseqnum,chainID))
    for r in R:
        if r.chainID==chainID and r.resseqnum==resseqnum and r.insertion==insertion:
            for a in r.atoms:
                if a.name==atname:
#                    print('returning',a)
                    return a
#            print('Error: no atom named {} found in {}{} of chain {}'.format(at.name,r.resseqnum,r.name,r.chainID))
#    print('Error: resid {} not found in chain {}'.format(resseqnum,chainID))
    return ''

if __name__=='__main__':
    r1=Residue()
    r1.resseqnum=381
    r1.insertion=' '
    r2=Residue()
    r2.resseqnum=381
    r2.insertion='A'
    print(r2<r1)