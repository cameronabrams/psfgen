from residue import _PDBResName123_

class Deletion:
    def __init__(self,commandlinerecord=''):
        if len(commandlinerecord)>0:
            self.commandlinerecord=commandlinerecord
            if not commandlinerecord[0].isalpha() or commandlinerecord[1]!='_':
                print('Poorly formed deltion spec: {}'.format(commandlinerecord))
                self.chainID='*'
                self.orig='***'
                self.new='***'
                self.orig_1='*'
                self.new_1='*'
                self.resseqnum=-999
            else:
                # C_Nrrr
                self.chainID=commandlinerecord[0]
                self.orig_1=commandlinerecord[2]
                self.orig=_PDBResName123_[self.orig_1]
                self.resseqnum=int(commandlinerecord[3:])
    def __str__(self):
        return '{}_{}{}'.format(self.chainID,self.orig_1,self.resseqnum)
    def replicate(self,newchainID=''):
        return Deletion(commandlinerecord=self.deletionStr(newChainID=newchainID))
    def deletionStr(self,newChainID=''):
        return '{}_{}{}'.format(self.chainID if newChainID=='' else newChainID,self.orig_1,self.resseqnum)


