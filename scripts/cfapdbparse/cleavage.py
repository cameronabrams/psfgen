class Cleavage:
    def __init__(self,commandlinerecord):
        self.commandlinerecord=commandlinerecord
        if not commandlinerecord[0].isalpha():
            print('Poorly formed cleavage spec: {}'.format(commandlinerecord))
            self.parent_chainID='*'
            self.parent_Cterm_resseqnum=-999
        else:
            self.parent_chainID=commandlinerecord[0]
            self.parent_Cterm_resseqnum=int(commandlinerecord[1:])
    def Clone(self,parent_chain=''):
        if len(parent_chain)>0:
            newCleavage=Cleavage(self.commandlinerecord)
            newCleavage.parent_chainID=parent_chain
            return newCleavage
    def __str__(self):
        return '{}{}-x-'.format(self.parent_chainID,self.parent_Cterm_resseqnum)


