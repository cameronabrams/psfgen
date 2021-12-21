from residue import _PDBResName123_

class Mutation:
    def __init__(self,commandlinerecord='',seqadv=''):
        if len(commandlinerecord)>0:
            self.commandlinerecord=commandlinerecord
            ''' if no chainID in the first character or the _ separator is not the second character '''
            if not commandlinerecord[0].isalpha() or commandlinerecord[1]!='_':
                print('Poorly formed mutation spec: {}'.format(commandlinerecord))
                self.chainID='*'
                self.orig='***'
                self.new='***'
                self.orig_1='*'
                self.new_1='*'
                self.resseqnumi='0 '
                self.resseqnum=0
                self.insertion=' '
            else:
                # (chainID:1s)_{orig:1s}{resseqnumi:s}{new:1s}
                self.chainID=commandlinerecord[0]
                self.orig_1=commandlinerecord[2] # single-letter AA designation
                self.orig=_PDBResName123_[self.orig_1]
                self.new_1=commandlinerecord[-1] # single-letter AA designation
                self.new=_PDBResName123_[self.new_1]
                self.resseqnumi=commandlinerecord[3:-1] # resnum+insertion
                if self.resseqnumi[-1].isalpha():
                    self.insertion=self.resseqnumi[-1]
                    self.resseqnum=int(self.resseqnumi[:-1])
                else:
                    self.insertion=' '
                    self.resseqnum=int(self.resseqnumi)
        elif seqadv!='':
            self.commandlinerecord=seqadv.pdbrecord
            self.chainID=seqadv.chainID
            self.orig=seqadv.resName
            self.new=seqadv.dbRes
            self.resseqnum=seqadv.seqNum
            self.insertion=seqadv.iCode
    def __str__(self):
        return f'{self.commandlinerecord} => chain {self.chainID} resseqnum {self.resseqnum} (insertion [{self.insertion}]) orig_resname {self.orig} mutant_resname {self.new}'
    def Clone(self,chain=''):
        if len(chain)>0:
            newMutation=Mutation(commandlinerecord=self.commandlinerecord)
            newMutation.chainID=chain
            newMutation.commandlinerecord=newMutation.mutationStr()
            return newMutation
    def mutationStr(self):
        return f'{self.chainID}_{self.orig_1}{self.resseqnumi}{self.new_1}'
    def psfgen_segment_str(self):
        return '   mutate {} {}\n'.format(self.resseqnumi,self.new) if self.chainID!='*' else ''

if __name__=='__main__':
    str1='C_L981AF'
    m1=Mutation(commandlinerecord=str1)
    print(str(m1))
    m2=m1.Clone(chain='V')
    print(str(m2))