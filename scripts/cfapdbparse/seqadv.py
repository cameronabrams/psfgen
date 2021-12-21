class Seqadv:
    def __init__(self,pdbrecord):
        if len(pdbrecord)>0:
            self.pdbrecord=pdbrecord
            self.record_name=pdbrecord[0:6]
            self.idCode=pdbrecord[7:11]
            self.rawresName=pdbrecord[12:15]
            self.resName=self.rawresName.strip()
            self.chainID=pdbrecord[16:17]
            self.seqNum=int(pdbrecord[18:22])
            self.iCode=pdbrecord[22:23]
            self.database=pdbrecord[24:28]
            self.dbAccession=pdbrecord[29:38]
            self.dbRes=pdbrecord[39:42].strip()
            self.dbSeqraw=pdbrecord[43:48]
            if self.dbSeqraw.strip().isdigit()==True:
                self.dbSeq=int(self.dbSeqraw)
            else:
                self.dbSeq=''
            self.conflict=pdbrecord[49:70].strip()
    def pdb_line(self):
        return f'{self.record_name:6s} {self.idCode:3s} {self.rawresName:>3s} {self.chainID:1s} {self.seqNum:>4d}{self.iCode:1s} {self.database:>4s} {self.dbAccession:9s} {self.dbRes:3s} {self.dbSeqraw:5s} {self.conflict:21s}'
    def Clone(self,chain=''):
        if len(chain)==1:
            newSeqadv=Seqadv(pdbrecord=self.pdb_line())
            newSeqadv.chainID=chain
            newSeqadv.pdbrecord=newSeqadv.pdb_line()
            return newSeqadv
    def printshort(self):
        return '{}-{}{}{}'.format(self.chainID,self.resName,self.seqNum,self.dbRes)
    def __str__(self):
        retstr='{}\n'+\
               ' idCode      {}\n'+\
               ' resName     {}\n'+\
               ' chainID     {}\n'+\
               ' seqNum      {}\n'+\
               ' iCode       {}\n'+\
               ' database    {}\n'+\
               ' dbAccession {}\n'+\
               ' dbRes       {}\n'+\
               ' dbSeq       {}\n'+\
               ' conflict    {}\n'
        return retstr.format(self.record_name,
                             self.idCode,
                             self.resName,
                             self.chainID,
                             self.seqNum,
                             self.iCode,
                             self.database,
                             self.dbAccession,
                             self.dbRes,
                             self.dbSeq,
                             self.conflict)

if __name__=='__main__':
    pr1='SEQADV 4ZMJ ASN G  332  UNP  Q2N0S6    THR   330 ENGINEERED MUTATION            '       
    pr2='SEQADV 4ZMJ CYS G  501  UNP  Q2N0S6    ALA   498 ENGINEERED MUTATION            '
    pr3='SEQADV 4ZMJ ARG G  508  UNP  Q2N0S6              EXPRESSION TAG                 '
    pr4='SEQADV 4ZMJ ARG G  509  UNP  Q2N0S6              EXPRESSION TAG                 '
    s1=Seqadv(pdbrecord=pr1)
    s2=Seqadv(pdbrecord=pr2)
    s3=Seqadv(pdbrecord=pr3)
    s4=Seqadv(pdbrecord=pr4)

    print(str(s1))
    print(pr1)
    print(s1.pdb_line())

    s5=s1.Copy(chain='F')
    print(str(s5))
