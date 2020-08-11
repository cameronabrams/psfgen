class Seqadv:
    def __init__(self,pdbrecord):
        if len(pdbrecord)>0:
            self.pdbrecord=pdbrecord
            self.record_name=pdbrecord[0:6]
            self.idCode=pdbrecord[7:11]
            self.resName=pdbrecord[12:15].strip()
            self.chainID=pdbrecord[16:17]
            self.seqNum=int(pdbrecord[18:22])
            self.iCode=pdbrecord[22:23]
            self.database=pdbrecord[24:28]
            self.dbAccession=pdbrecord[29:38]
            self.dbRes=pdbrecord[39:42].strip()
            if pdbrecord[43:48].strip().isdigit()==True:
                self.dbSeq=int(pdbrecord[43:48])
            else:
                self.dbSeq=''
            self.conflict=pdbrecord[49:70].strip()
    def pdbline(self):
            print(self.pdbrecord)
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

