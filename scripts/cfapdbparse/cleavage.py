class Cleavage:
    def __init__(self,commandlinerecord):
        self.commandlinerecord=commandlinerecord
        if not commandlinerecord[0].isalpha() or not commandlinerecord[-1].isalpha():
            print('Poorly formed cleavage spec: {}'.format(commandlinerecord))
            self.parent_chainID='*'
            self.parent_Cterm_resseqnum=-999
            self.daughter_chainID='*'
        else:
            self.parent_chainID=commandlinerecord[0]
            self.parent_Cterm_resseqnum=int(commandlinerecord[1:-1])
            self.daughter_chainID=commandlinerecord[-1]
    def __str__(self):
        return '{}{}-x-{}'.format(self.parent_chainID,self.parent_Cterm_resseqnum,self.daughter_chainID)


