from datetime import datetime
class RevDat:
    def __init__(self,pdbrecord):
        self.pdbrecord=pdbrecord
        self.datecode=datetime.strptime(pdbrecord[13:22],'%d-%b-%y')
        self.pdbcode=pdbrecord[23:27]
        self.isnotfirs=int(pdbrecord[31:32])
        self.tokens=self.get_tokens(pdbrecord)
    def __str__(self):
        return self.datecode.strftime('%B %d, %Y')
    def show(self):
        retstr=self.datecode.strftime('%B %d, %Y: ')
        for t in self.tokens:
            retstr+=t+' '
        return retstr
    def addTokens(self,pdbrecord):
        ctok=pdbrecord[11:12]
        if ctok.isdigit():
            self.tokens.extend(self.get_tokens(pdbrecord))
    def get_tokens(self,pdbrecord):
        tks=[]
        for t in pdbrecord[39:].strip().split(' '):
            tks.append(t.strip())
        return tks

class FmtDat:
    def __init__(self,pdbrecord):
        self.PDBVersion=float(pdbrecord[40:44])
        self.PDBVersionDate=datetime.strptime(pdbrecord[46:55],'%d-%b-%y')
    def __str__(self):
        return 'Format version {:.2f}, {}'.format(self.PDBVersion,self.PDBVersionDate.strftime('%B %d, %Y'))

