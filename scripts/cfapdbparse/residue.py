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
    def __init__(self,a=-1,m=0,molid='*'):
        ''' initializing with an atom '''
        if a!=-1:
            self.resseqnum=a.resseqnum
            self.insertion=a.insertion
            self.name=a.resname
            self.chainID=a.chainID
            self.source_chainID=a.chainID
            self.atoms=[a]
            self.up=[]
            self.down=[]
        else:
            ''' initializing with a missing residue '''
            if m!=0:
                self.resseqnum=m.resseqnum
                self.name=m.resname
                self.insertion=m.insertion
                self.chainID=m.chainID
                self.source_chainID=m.chainID
                self.biomt='*'
                self.atoms=[]
                self.up=[]
                self.down=[]
            else:
                print('ERROR: bad residue construction')

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
    def set_connections(self):
        pass
    def __str__(self):
        return '{}-{}{}'.format(self.chainID,self.name,self.resseqnum)
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
    for r in R:
        if r.chainID==chainID and r.resseqnum==resseqnum and r.insertion==insertion:
            return r
    return '' 

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

