from residue import _PDBResName123_

class Mutation:
    def __init__(self,commandlinerecord='',seqadv=''):
        if len(commandlinerecord)>0:
            self.commandlinerecord=commandlinerecord
            if not commandlinerecord[0].isalpha() or commandlinerecord[1]!='_':
                print('Poorly formed mutation spec: {}'.format(commandlinerecord))
                self.chainID='*'
                self.orig='***'
                self.new='***'
                self.orig_1='*'
                self.new_1='*'
                self.resseqnum=-999
            else:
                self.chainID=commandlinerecord[0]
                self.orig_1=commandlinerecord[2]
                self.orig=_PDBResName123_[self.orig_1]
                self.new_1=commandlinerecord[-1]
                self.new=_PDBResName123_[self.new_1]
                self.resseqnum=int(commandlinerecord[3:-1])
        elif seqadv!='':
            self.commandlinerecord=seqadv.pdbrecord
            self.chainID=seqadv.chainID
            self.orig=seqadv.resName
            self.new=seqadv.dbRes
            self.resseqnum=seqadv.seqNum
    def __str__(self):
        return '{}_{}{}{}'.format(self.chainID,self.orig_1,self.resseqnum,self.new_1)
    def replicate(self,newchainID=''):
        return Mutation(commandlinerecord=self.mutationStr(newChainID=newchainID))
    def mutationStr(self,newChainID=''):
        return '{}_{}{}{}'.format(self.chainID if newChainID=='' else newChainID,self.orig_1,self.resseqnum,self.new_1)
    def psfgen_segment_str(self):
        return '   mutate {} {}\n'.format(self.resseqnum,self.new) if self.chainID!='*' else ''

