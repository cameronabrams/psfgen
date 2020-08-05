class Missing:
    def __init__(self,pdbrecord):
        self.pdbrecord=pdbrecord
        self.record_name=pdbrecord[0:6]
        self.code=int(pdbrecord[7:10])
        self.model=pdbrecord[13:14].strip()
        self.resname=pdbrecord[15:18].strip()
        self.chainID=pdbrecord[19:20]
        self.resseqnum=int(pdbrecord[21:26])
    def __str__(self):
        retstr='MISSING\n'+\
               '   model     {:s}\n'+\
               '   resname   {:s}\n'+\
               '   chainID   {:s}\n'+\
               '   resseqnum {:d}\n'
        return retstr.format(self.model,self.resname,self.chainID,self.resseqnum)
    def psfgen_residueline(self):
        return '     residue {} {} {}'.format(self.resname,self.resseqnum,self.chainID)

