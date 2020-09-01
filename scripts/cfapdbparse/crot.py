class Crot:
    ''' user-defined rotation of phi/psi/chiX bonds of
        specific residue to reorient all atoms C-terminal
        to that residue up to a designated end residue 
    '''
    def __init__(self,record):
        self.record=record
        dat=record.split(',')
        self.angle=dat[0].upper()
        if self.angle=='PHI' or self.angle=='PSI' or self.angle=='OMEGA':
            # this is a backbone torsion, so we need both an owner
            # residue and a residue indicating the end of the 
            # set of residues that will be reoriented by the
            # rotation
            self.chainID=dat[1]
            self.resseqnum1=int(dat[2])
            self.resseqnum2=int(dat[3])
            self.degrees=float(dat[4])
        elif self.angle=='CHI1' or self.angle=='CHI2':
            self.chainID=dat[1]
            self.resseqnum1=int(dat[2])
            self.resseqnum2=-1
            self.degrees=float(dat[3])
        else:
            print('Error: unable to parse Crot argument {}'.format(record))

    def __str__(self):
        return self.record

    def psfgen_str(self,molid='top'):
        retstr=''
        if self.resseqnum2!=-1:  # this is a backbone bond
            retstr+='set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]\n'.format(molid,self.chainID,self.resseqnum1)
            retstr+='set r2 [[atomselect {} "chain {} and resid {} and name CA"] get residue]\n'.format(molid,self.chainID,self.resseqnum2)
            retstr+='Crot_{} $r1 $r2 {} {} {}\n'.format(self.angle.lower(),self.chainID,molid,self.degrees)
        elif self.resseqnum2==-1:  # this is a side-chain bond
            retstr+='set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]\n'.format(molid,self.chainID,self.resseqnum1)
            retstr+='SCrot_{} $r1 {} {} {}\n'.format(self.angle.lower(),self.chainID,molid,self.degrees)
        else:
            print('Warning: buggy crot {}'.format(self.record))
        return retstr

        