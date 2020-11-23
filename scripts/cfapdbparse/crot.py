class Crot:
    ''' user-defined rotation and translations
    '''
    def __init__(self,record):
        self.record=record
        dat=record.split(',')
        self.angle=dat[0].upper()
#        print(dat)
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
        elif self.angle=='GLYCAN':
            self.segname=dat[1]
            self.resseqnum1=int(dat[2])
            self.atom1=dat[3]
            self.resseqnum2=int(dat[4])
            self.atom2=dat[5]
            self.degrees=float(dat[6])
        elif self.angle=='LINK':
            self.segname1=dat[1]
            self.resseqnum1=int(dat[2])
            self.atom1=dat[3]
            self.segname2=dat[4]
            self.resseqnum2=int(dat[5])
            self.atom2=dat[6]
            self.degrees=float(dat[7])
        elif self.angle=='ANGLEIJK':
            self.segnamei=dat[1]
            self.resseqnumi=int(dat[2])
            self.atomi=dat[3]
            self.segnamejk=dat[4]
            self.resseqnumj=int(dat[5])
            self.atomj=dat[6]
            self.resseqnumk=int(dat[7])
            self.atomk=dat[8]
            self.degrees=float(dat[9])
        else:
            print('Error: unable to parse Crot argument {}'.format(record))
    def replicate(self,newc='',newsegname=''):
        newcrot=Crot(self.record)
        if newc!='':
            newcrot.chainID=newc
        if newsegname!='':
            newcrot.segname=newsegname
        return newcrot
    def __str__(self):
        return self.record

    def psfgen_str(self,molid='top'):
        retstr=''
        if self.angle=='PHI' or self.angle=='PSI' or self.angle=='OMEGA':  # this is a backbone bond
            retstr+='set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]\n'.format(molid,self.chainID,self.resseqnum1)
            retstr+='set r2 [[atomselect {} "chain {} and resid {} and name CA"] get residue]\n'.format(molid,self.chainID,self.resseqnum2)
            retstr+='Crot_{} $r1 $r2 {} {} {}\n'.format(self.angle.lower(),self.chainID,molid,self.degrees)
        elif self.angle=='CHI1' or self.angle=='CHI2':  # this is a side-chain bond
            retstr+='set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]\n'.format(molid,self.chainID,self.resseqnum1)
            retstr+='SCrot_{} $r1 {} {} {}\n'.format(self.angle.lower(),self.chainID,molid,self.degrees)
        elif self.angle=='GLYCAN':  # intra-glycan rotation
            retstr+='set sel [atomselect {} "segname {}"]\n'.format(molid,self.segname1)
            retstr+='set i [[atomselect {} "segname {} and resid {} and name {}"] get index]\n'.format(molid,self.segname1,self.resseqnum1,self.atom1)
            retstr+='set j [[atomselect {} "segname {} and resid {} and name {}"] get index]\n'.format(molid,self.segname2,self.resseqnum2,self.atom2)
            retstr+='genbondrot {} $sel $i $j {}\n'.format(molid,self.degrees)
        elif self.angle=='LINK': # ASN-GLYcan rotation
            retstr+='set sel [atomselect {} "segname {} {}"]\n'.format(molid,self.segname1,self.segname2)
            retstr+='set i [[atomselect {} "segname {} and resid {} and name {}"] get index]\n'.format(molid,self.segname1,self.resseqnum1,self.atom1)
            retstr+='set j [[atomselect {} "segname {} and resid {} and name {}"] get index]\n'.format(molid,self.segname2,self.resseqnum2,self.atom2)
            retstr+='genbondrot {} $sel $i $j {}\n'.format(molid,self.degrees)
        elif self.angle=='ANGLEIJK':
            retstr+='set rotsel [atomselect {} "segname {}"]\n'.format(molid,self.segnamejk)
            retstr+='set ri [lindex [atomselect {} "segname {} and resid {} and name {}] get {x y z}] 0]\n'.format(molid,self.segnamei,self.reseqnumi,self.atomi)
            retstr+='set rj [lindex [atomselect {} "segname {} and resid {} and name {}] get {x y z}] 0]\n'.format(molid,self.segnamej,self.reseqnumj,self.atomj)
            retstr+='set rk [lindex [atomselect {} "segname {} and resid {} and name {}] get {x y z}] 0]\n'.format(molid,self.segnamek,self.reseqnumk,self.atomk)
            retstr+='set rij [vecsub $ri $rj]\n'
            retstr+='set rjk [vecsub $rj $rk]\n'
            retstr+='set cijk [veccross $rij $rjk]\n'
            retstr+='$rotsel move [trans center $rj origin $rj axis $cijk {} degrees]\n'.format(self.degrees)
        else:
            print('Warning: buggy crot {}'.format(self.record))
        return retstr

        
