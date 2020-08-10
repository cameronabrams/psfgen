from residue import _PDBResName123_

class Mutation:
    def __init__(self,commandlinerecord):
        self.commandlinerecord=commandlinerecord
        if not commandlinerecord[0].isalpha() or commandlinerecord[1]!='_':
            print('Poorly formed mutation spec: {}'.format(commandlinerecord))
            self.chainID='*'
            self.orig='*'
            self.new='*'
            self.resseqnum=-999
        else:
            self.chainID=commandlinerecord[0]
            self.orig=_PDBResName123_[commandlinerecord[2]]
            self.new=_PDBResName123_[commandlinerecord[-1]]
            self.resseqnum=int(commandlinerecord[3:-1])
    def __str__(self):
        return '{}-{}{}{}'.format(self.chainID,self.orig,self.resseqnum,self.new)
    def psfgen_segment_str(self):
        return '   mutate {} {}\n'.format(self.resseqnum,self.new) if self.chainID!='*' else ''
