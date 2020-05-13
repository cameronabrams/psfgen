_PDBAtomNameDict_={'CL':'CLA'}
class Atom:
    def __init__(self,line=[]):
        if len(line)>0:
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
#31 - 38       Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
            x=float(line[30:38])
#39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
            y=float(line[38:46])
#47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
            z=float(line[46:54])
#55 - 60        Real(6.2)     occupancy    Occupancy.
            occ=float(line[54:60])
#61 - 66       Real(6.2)     tempFactor   Temperature  factor.
            beta=float(line[60:66])
#77 - 78        LString(2)    element      Element symbol, right-justified.
            elem=line[76:78]
#79 - 80        LString(2)    charge       Charge  on the atom.
            charge=line[78:80]
            self.record_name=record_name.strip()
            self.serial=serial
            self.name=name.strip()
            self.altloc=altloc.strip()
            self.resname=resname.strip()
            self.chainID=chainID.strip()
            self.resseqnum=resseqnum
            self.insertion=insertion.strip()
            self.x=x
            self.y=y
            self.z=z
            self.occ=occ
            self.beta=beta
            self.elem=elem.strip()
            self.charge=charge.strip()
            self.empty=False
        else:
            self.empty=True 
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


