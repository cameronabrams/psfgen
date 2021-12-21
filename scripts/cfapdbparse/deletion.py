from residue import _PDBResName123_

class Deletion:
    def __init__(self,commandlinerecord=''):
        if len(commandlinerecord)>0:
            self.commandlinerecord=commandlinerecord
            if not commandlinerecord[0].isalpha() or commandlinerecord[1]!='_':
                print('Poorly formed deltion spec: {}'.format(commandlinerecord))
                self.chainID='*'
                self.resseqnumi='0 '
                self.resseqnum=0
                self.insertion=' '
            else:
                # C_Nrrr
                self.chainID=commandlinerecord[0]
                self.resseqnumi=commandlinerecord[3:]
                if self.resseqnumi[-1].isalpha():
                    self.insertion=self.resseqnumi[-1]
                    self.resseqnum=int(self.resseqnumi[:-1])
                else:
                    self.resseqnum=int(self.resseqnumi)
                    self.insertion=' '
    def __lt__(self,other):
        if self.resseqnum<other.resseqnum:
            return True
        elif self.resseqnum==other.resseqnum:
            if self.insertion==None and other.insertion==None:
                return False
            elif (self.insertion=='' or self.insertion==' ' or self.insertion==None) and other.insertion.isalpha():
                return True
            elif self.insertion.isalpha() and other.insertion.isalpha():
                return ord(self.insertion)<ord(other.insertion)
            else:
                return False
    def __str__(self):
        return f'{self.commandlinerecord} => chain {self.chainID} resseqnum {self.resseqnum} (insertion [{self.insertion}])'
    def Clone(self,chain=''):
        if len(chain)>1:
            newDeletion=Deletion(commandlinerecord=self.commandlinerecord)
            newDeletion.chainID=chain
            return newDeletion
#    def replicate(self,newchainID=''):
#        return Deletion(commandlinerecord=self.deletionStr(newChainID=newchainID))
    def deletionStr(self,newChainID=''):
        return '{}_{}{}'.format(self.chainID if newChainID=='' else newChainID,self.orig_1,self.resseqnum)

if __name__=='__main__':
    from residue import Residue
    r1=Residue()
    r1.resseqnum=381
    r1.insertion=' '
    r2=Residue()
    r2.resseqnum=381
    r2.insertion='A'
    print(r1<r2)
    d1=Deletion(commandlinerecord='C_381')
    print(r1<d1)
    print(r2<d1)
