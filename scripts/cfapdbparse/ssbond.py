class SSBond:
    def __init__(self,pdbrecord=[]):
# 1 -  6        Record name    "SSBOND"
        record_name=pdbrecord[0:6]
# 8 - 10        Integer        serNum           Serial number.
        serial_number=int(pdbrecord[7:10])
# 12 - 14        LString(3)     "CYS"            Residue name.
        resname1=pdbrecord[11:14]
# 16             Character      chainID1         Chain identifier.
        chainID1=pdbrecord[15:16]
# 18 - 21        Integer        seqNum1          Residue sequence number.
        resseqnum1=int(pdbrecord[17:21])
# 22             AChar          icode1           Insertion code.
        icode1=pdbrecord[21:22]
# 26 - 28        LString(3)     "CYS"            Residue name.
        resname2=pdbrecord[25:28]
# 30             Character      chainID2         Chain identifier.
        chainID2=pdbrecord[29:30]
# 32 - 35        Integer        seqNum2          Residue sequence number.
        resseqnum2=int(pdbrecord[31:35])
# 36             AChar          icode2           Insertion code.
        icode2=pdbrecord[35:36]
# 60 - 65        SymOP          sym1             Symmetry operator for residue 1.
        sym1=pdbrecord[59:65]
# 67 - 72        SymOP          sym2             Symmetry operator for residue 2.
        sym2=pdbrecord[66:72]
# 74 – 78        Real(5.2)      Length           Disulfide bond distance
        length=float(pdbrecord[73:78])

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
       return 'patch DISU {}:{} {}:{}\n'.format(self.chainID1,self.resseqnum1,self.chainID2,self.resseqnum2)


