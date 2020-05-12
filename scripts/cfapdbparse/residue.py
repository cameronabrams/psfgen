_pdb_ions_=['LIT','SOD','MG','POT','CAL','RUB','CES','BAR','ZN','CAD','CL']
_pdb_glycans_=['BMA','FUC','GAL','MAN','NAG','SIA']

_PDBResName123_={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY',
               'H':'HSE','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER',
               'T':'THR','W':'TRP','Y':'TYR','V':'VAL'}
_PDBResNameDict_={'HIS':'HSE','ZN':'ZN2','HOH':'TIP3','CL':'CLA','NAG':'BGNA','MAN':'AMAN','BMA':'BMAN','FUC':'AFUC','GAL':'BGAL','SIA':'ANE5AC'}

class Residue:
    def __init__(self,a=-1,m=0):
        if a!=-1:
            self.resseqnum=a.resseqnum
            self.name=a.resname
            self.chainID=a.chainID
            self.source_chainID=a.chainID
            self.atoms=[a]
        else:
            if m!=0:
                self.resseqnum=m.resseqnum
                self.name=m.resname
                self.chainID=m.chainID
                self.source_chainID=m.chainID
                self.atoms=[]
            else:
                fp.write('ERROR: bad residue construction')

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
    def __str__(self):
        if len(self.atoms)==0:
            atstr='MISSING'
        else:
            atser=[]
            for a in self.atoms:
                atser.append(a.serial)
            atstr='{:d} - {:d}'.format(min(atser),max(atser))
        return 'RESIDUE {} {} {:d} {}'.format(self.chainID,self.name,self.resseqnum,atstr)


