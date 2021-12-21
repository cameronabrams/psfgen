class Missing:
    def __init__(self,pdbrecord=None,cifdict=None):
        if pdbrecord!=None and cifdict!=None:
            print("Error: <Missing> __init__ called with both a pdbrecord and cifdict.\nUsing pdbrecord.")
        if pdbrecord!=None:
            self.pdbrecord=pdbrecord
            self.record_name=pdbrecord[0:6]
            self.code=int(pdbrecord[7:10])
            self.rawmodel=pdbrecord[13:14]
            self.model=self.rawmodel.strip()
            self.rawresname=pdbrecord[15:18]
            self.resname=self.rawresname.strip()
            self.chainID=pdbrecord[19:20]
            self.resseqnum=int(pdbrecord[21:26])
            self.insertion=pdbrecord[26:27]
        elif cifdict!=None:
            self.cifdict=cifdict
            self.model=int(cifdict['pdb_model_num'])
            self.resname=cifdict['auth_comp_id']
            self.chainID=cifdict['auth_asym_id']
            self.resseqnum=int(cifdict['auth_seq_id'])
            ic=cifdict['pdb_ins_code']
            self.insertion=' ' if ic=='?' else ic
    def pdb_line(self):
        return '{:6s}{:>4d}   {:1s} {:3s} {:1s} {:>5d}{:1s}'.format(self.record_name,
        self.code,self.rawmodel,self.rawresname,self.chainID,self.resseqnum,self.insertion)
    def Clone(self,chain=''):
        if len(chain)==1:
            newMissing=Missing(pdbrecord=self.pdb_line())
            newMissing.chainID=chain
            newMissing.pdbrecord=newMissing.pdb_line()
            return newMissing
    def __str__(self,verbose=False):
        if verbose:
            retstr='MISSING\n'+\
                   '   model     {:s}\n'+\
                   '   resname   {:s}\n'+\
                   '   chainID   {:s}\n'+\
                   '   resseqnum {:d}\n'+\
                   '   insertion {:s}\n'
        else:
            retstr='{}{}_{}{}{}'
        return retstr.format(self.model,self.resname,self.chainID,self.resseqnum,self.insertion)
    def psfgen_residueline(self):
        return '     residue {} {}{} {}'.format(self.resname,self.resseqnum,self.insertion,self.chainID)

if __name__=='__main__':
    pr='REMARK 465     GLU G   185A  '
    m1=Missing(pdbrecord=pr)
    print(str(m1))
    print(pr)
    print(m1.pdb_line())
    m2=m1.Clone(chain='Z')
    print(str(m2))