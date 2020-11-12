class Missing:
    def __init__(self,pdbrecord=None,cifdict=None):
        if pdbrecord!=None and cifdict!=None:
            print("Error: Missing __init__ called with both a pdbrecord and cifdict.\nUsing pdbrecord.")
        if pdbrecord!=None:
            self.pdbrecord=pdbrecord
            self.record_name=pdbrecord[0:6]
            self.code=int(pdbrecord[7:10])
            self.model=pdbrecord[13:14].strip()
            self.resname=pdbrecord[15:18].strip()
            self.chainID=pdbrecord[19:20]
            self.resseqnum=int(pdbrecord[21:26])
            self.insertion=pdbrecord[26:27]
        elif cifdict!=None:
            self.cifdict=cifdict
            self.model=int(cifdict['pdb_model_num'])
            self.resname=cifdict['auth_comp_id']
            self.chainID=cifdcit['auth_asym_id']
            self.resseqnum=int(cifdict['auth_seq_id'])
            ic=cifdict['pdb_ins_code']
            self.insertion=' ' if ic=='?' else ic
    def __str__(self):
        retstr='MISSING\n'+\
               '   model     {:s}\n'+\
               '   resname   {:s}\n'+\
               '   chainID   {:s}\n'+\
               '   resseqnum {:d}\n'+\
               '   insertion {:s}\n'
        return retstr.format(self.model,self.resname,self.chainID,self.resseqnum,self.insertion)
    def psfgen_residueline(self):
        return '     residue {} {}{} {}'.format(self.resname,self.resseqnum,self.insertion,self.chainID)

