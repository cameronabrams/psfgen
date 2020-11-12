class SSBond:
    ''' 
    Creates SSBond instance from one of two types of input strings:
    1. a PDB SSBOND record
    2. a command-line argument of the form X_###-Y_### where X and Y are chain IDs and
    the ###'s are resids
    Can also accept a CIFDict
    '''
    def __init__(self,record=None,cifdict=None):
        if record!=None and cifdict!=None:
            print('Error: SSBond __init__ called with both a record and a cifdict.\nUsing the record.')
        if record!=None:
            if 'SSBOND' in record:
                self.ssbond_from_pdbrecord(record)
            else: # assume this is a command-line argument
                self.ssbond_from_commandline(record)
        elif cifdict!=None:
            self.ssbond_from_cifdict(cifdict)
    def ssbond_from_cifdict(self,d):
        self.record_name='SSBOND'
        self.serial_number=int(d['id'].strip('disulf'))
        self.resname1='CYS'
        self.resname2='CYS'
        self.chainID1=d['ptnr1_auth_asym_id']
        self.chainID2=d['ptnr2_auth_asym_id']
        self.resseqnum1=int(d['ptnr1_auth_seq_id'])
        self.resseqnum2=int(d['ptnr2_auth_seq_id'])
        self.icode1=d['pdbx_ptnr1_pdb_ins_code']
        self.icode2=d['pdbx_ptnr2_pdb_ins_code']
        self.sym1=d['ptnr1_symmetry']
        self.sym2=d['ptnr2_symmetry']
        self.length=float(d['pdbx_dist_value'])
        self.pdbrecord=self.pdb_line()
    def ssbond_from_commandline(self,record):
        self.record_name='SSBOND'
        self.serial_number=0
        self.resname1='CYS'
        self.resname2='CYS'
        s1=record.split('-')
        r1=s1[0].split('_')
        r2=s1[1].split('_')
        self.chainID1=r1[0]
        self.chainID2=r2[0]
        self.resseqnum1=r1[1]
        self.resseqnum2=r2[1]
        self.icode1=''
        self.icode2=''
        self.length=0.0
        self.sym1=''
        self.sym2=''
        self.pdbrecord=self.pdb_line()
    def ssbond_from_pdbrecord(self,pdbrecord):
        self.pdbrecord=pdbrecord
# 1 -  6        Record name    "SSBOND"
        self.record_name=pdbrecord[0:6].strip()
# 8 - 10        Integer        serNum           Serial number.
        self.serial_number=int(pdbrecord[7:10])
# 12 - 14        LString(3)     "CYS"            Residue name.
        self.resname1=pdbrecord[11:14].strip()
# 16             Character      chainID1         Chain identifier.
        self.chainID1=pdbrecord[15:16].strip()
# 18 - 21        Integer        seqNum1          Residue sequence number.
        self.resseqnum1=int(pdbrecord[17:21])
# 22             AChar          icode1           Insertion code.
        self.icode1=pdbrecord[21:22].strip()
# 26 - 28        LString(3)     "CYS"            Residue name.
        self.resname2=pdbrecord[25:28].strip()
# 30             Character      chainID2         Chain identifier.
        self.chainID2=pdbrecord[29:30].strip()
# 32 - 35        Integer        seqNum2          Residue sequence number.
        self.resseqnum2=int(pdbrecord[31:35])
# 36             AChar          icode2           Insertion code.
        self.icode2=pdbrecord[35:36].strip()
# 60 - 65        SymOP          sym1             Symmetry operator for residue 1.
        self.sym1=pdbrecord[59:65].strip()
# 67 - 72        SymOP          sym2             Symmetry operator for residue 2.
        self.sym2=pdbrecord[66:72].strip()
# 74 – 78        Real(5.2)      Length           Disulfide bond distance
        self.length=float(pdbrecord[73:78])
    def pdb_line(self):
        pdbline='{:6s}'.format(self.record_name)+\
                '{:4d}'.format(self.serial_number)+\
                '{:>4s}'.format(self.resname1)+\
                '{:>2s}'.format(self.chainID1)+\
                '{:5d}'.format(self.resseqnum1)+\
                '{:1s}'.format(self.icode1)+\
                '{:>6s}'.format(self.resname2)+\
                '{:>2s}'.format(self.chainID2)+\
                '{:5d}'.format(self.resseqnum2)+\
                '{:1s}'.format(self.icode2)+\
                ' '*23+\
                '{:>6s}'.format(self.sym1)+\
                '{:>7s}'.format(self.sym2)+\
                '{:6.2f}'.format(self.length)
        return pdbline
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


