_PDBAtomNameDict_={'CL':'CLA'}
class Atom:
    def __init__(self,pdbrecord=None,cifdict=None):
        if pdbrecord!=None and cifdict!=None:
            print('Error: Atom __init__() called with both a pdbrecord and a cifdict.\nUsing the pdbrecord.')
        self.pdbrecord=pdbrecord
        self.cifdict=cifdict

        if pdbrecord!=None:
# 1 -  6        Record name   "ATOM  "
            self.record_name=pdbrecord[0:6].strip()
# 7 - 11        Integer       serial       Atom  serial number.
            self.serial=int(pdbrecord[6:11])
#13 - 16        Atom          name         Atom name.
            self.name=pdbrecord[12:16].strip() # modified for 4-char atom names
#17             Character     altLoc       Alternate location indicator.
            self.altloc=pdbrecord[16:17]
#18 - 20        Residue name  resName      Residue name.
            self.resname=pdbrecord[17:21].strip() # modified for 4-char resnames
#22             Character     chainID      Chain identifier.
            self.chainID=pdbrecord[21:22]   # only one-character chain IDs are allowed
#23 - 26        Integer       resSeq       Residue sequence number.
            self.resseqnum=int(pdbrecord[22:26])
#27             AChar         iCode        Code for insertion of residues.
            self.insertion=pdbrecord[26:27]
#31 - 38       Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
            self.x=float(pdbrecord[30:38])
#39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
            self.y=float(pdbrecord[38:46])
#47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
            self.z=float(pdbrecord[46:54])
#55 - 60        Real(6.2)     occupancy    Occupancy.
            self.occ=float(pdbrecord[54:60])
#61 - 66       Real(6.2)     tempFactor   Temperature  factor.
            self.beta=float(pdbrecord[60:66])
#77 - 78        LString(2)    element      Element symbol, right-justified.
            self.elem=pdbrecord[76:78].strip()
#79 - 80        LString(2)    charge       Charge  on the atom.
            self.charge=pdbrecord[78:80].strip()
            self.segname=self.chainID
            self.empty=False
            self.link='None'
#            self.biomt=0
        elif cifdict!=None:
            self.serial=int(cifdict['id'])
            self.name=cifdict['auth_atom_id']
            al=cifdict['label_alt_id']
            self.altloc=' ' if al=='.' else al
            self.resname=cifdict['auth_comp_id']
            self.chainID=cifdict['auth_asym_id']
            self.resseqnum=int(cifdict['auth_seq_id'])
            ic=cifdict['pdbx_pdb_ins_code']
            self.insertion=' ' if ic=='?' else ic
            self.x=float(cifdict['cartn_x'])
            self.y=float(cifdict['cartn_y'])
            self.z=float(cifdict['cartn_z'])
            self.occ=float(cifdict['occupancy'])
            self.beta=float(cifdict['b_iso_or_equiv'])
            c=cifdict['pdbx_formal_charge']
            self.charge=0.0 if c=='?' else float(c)
            self.elem=cifdict['type_symbol']
            self.segname=self.chainID
            self.empty=False
            self.link='None'
        else:
            self.empty=True 
    def pdb_line(self):
        pdbline='{:<6s}'.format(self.record_name)+\
                '{:5d}'.format(self.serial)+' '+\
                '{:<4s}'.format(' '+self.name if len(self.name)<4 else self.name)+\
                '{:1s}'.format(self.altloc)+\
                '{:<4s}'.format(self.resname)+\
                '{:1s}'.format(self.chainID)+\
                '{:4d}'.format(self.resseqnum)+\
                '{:1s}'.format(self.insertion)+'   '+\
                '{:8.3f}'.format(self.x)+\
                '{:8.3f}'.format(self.y)+\
                '{:8.3f}'.format(self.z)+\
                '{:6.2f}'.format(self.occ)+\
                '{:6.2f}'.format(self.beta)+\
                10*' '+'{:>2s}'.format(self.elem)+'{:2s}'.format(self.charge)
        return pdbline
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


