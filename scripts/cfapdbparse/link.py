from segment import _seg_class_
class LinkSet:
    def __init__(self):
       self.L=[]
       self.C=[]
    def add_link(self,l):
       self.L.append(l)
    def cluster(self,Residues):
       for l in self.L:
           pass

class Link:
    def __init__(self,pdbrecord=[]):
# 1 -  6         Record name    "LINK  "
        if len(pdbrecord)>0:
            record_name=pdbrecord[0:6]
#13 - 16         Atom           name1           Atom name.
            name1=pdbrecord[12:16]
#17              Character      altLoc1         Alternate location indicator.
            altloc1=pdbrecord[16:17]
#18 - 20         Residue name   resName1        Residue  name.
            resname1=pdbrecord[17:20]
#    22              Character      chainID1        Chain identifier.
            chainID1=pdbrecord[21:22]
#23 - 26         Integer        resSeq1         Residue sequence number.
            resseq1=int(pdbrecord[22:26])
#27              AChar          iCode1          Insertion code.
            icode1=pdbrecord[26:27]
#43 - 46         Atom           name2           Atom name.
            name2=pdbrecord[42:46]
#47              Character      altLoc2         Alternate location indicator.
            altloc2=pdbrecord[46:47]
#48 - 50         Residue name   resName2        Residue name.
            resname2=pdbrecord[47:50]
#52              Character      chainID2        Chain identifier.
            chainID2=pdbrecord[51:52]
#53 - 56         Integer        resSeq2         Residue sequence number.
            resseq2=int(pdbrecord[52:56])
#57              AChar          iCode2          Insertion code.
            icode2=pdbrecord[56:57]  
#60 - 65         SymOP          sym1            Symmetry operator atom 1.
            sym1=pdbrecord[59:65]
#67 - 72         SymOP          sym2            Symmetry operator atom 2.
            sym2=pdbrecord[66:72]
#74 – 78         Real(5.2)      Length          Link distance
            link_distance=float(pdbrecord[73:78])
            self.record_name=record_name.strip()
            self.name1=name1.strip()
            self.altloc1=altloc1.strip()
            self.resname1=resname1.strip()
            self.chainID1=chainID1.strip()
            self.resseq1=resseq1
            self.icode1=icode1.strip()
            self.name2=name2.strip()
            self.altloc2=altloc2.strip()
            self.resname2=resname2.strip()
            self.chainID2=chainID2.strip()
            self.resseq2=resseq2
            self.icode2=icode2.strip()
            self.sym1=sym1.strip()
            self.sym2=sym2.strip()
            self.link_distance=link_distance
            self.empty=False
        else:
            self.empty=True

    def __str__(self):
        retstr='{}\n'+\
                '   name1         {:s}\n'+\
                '   altloc1       {:s}\n'+\
                '   resname1      {:s}\n'+\
                '   chainID1      {:s}\n'+\
                '   resseq1       {:d}\n'+\
                '   icode1        {:s}\n'+\
                '   name2         {:s}\n'+\
                '   altloc2       {:s}\n'+\
                '   resname2      {:s}\n'+\
                '   chainID2      {:s}\n'+\
                '   resseq2       {:d}\n'+\
                '   icode2        {:s}\n'+\
                '   sym1          {:s}\n'+\
                '   sym2          {:s}\n'+\
                '   link_distance {:.3f}'
        return retstr.format(self.record_name,self.name1,self.altloc1,self.resname1,self.chainID1,self.resseq1,self.icode1,self.name2,self.altloc2,self.resname2,self.chainID2,self.resseq2,self.icode2,self.sym1,self.sym2,self.link_distance)
    def psfgen_patchline(self):
        if self.resname1=='ASN' and _seg_class_[self.resname2]=='GLYCAN':
            return 'patch NGLB {}:{} {}S:{}\n'.format(self.chainID1,self.resseq1,self.chainID2,self.resseq2)
        else:
            # for a glycan-glycan patch, the C1 atom is always on the ji-residue
            if self.name2=='C1' and _seg_class_[self.resname1]=='GLYCAN':
                cmdj='[axeq {} 0 {} {} {}]'.format(self.resseq2,self.chainID2,self.name2,self.resseq1)
                cmdi='[axeq {} 0 {} {} {}]'.format(self.resseq1,self.chainID1,self.name1,-1)
                return 'patch 1{:1s}{}{} {}S:{} {}S:{}\n'.format(self.name1[1], cmdi,cmdj,self.chainID1,self.resseq1,self.chainID2,self.resseq2)
            elif self.name1=='C1' and _seg_class_[self.resname2]=='GLYCAN':
                cmdi='[axeq {} 0 {} {} {}]'.format(self.resseq2,self.chainID2,self.name2,self.resseq1)
                cmdj='[axeq {} 0 {} {} {}]'.format(self.resseq1,self.chainID1,self.name1,-1)           
                return 'patch 1{:1s}{}{} {}S:{} {}S:{}\n'.format(self.name2[1], cmdi,cmdj,self.chainID2,self.resseq2,self.chainID1,self.resseq1)
            elif self.name1=='O6' and self.name2=='C2':
                return 'patch SA26E {}S:{} {}S:{}\n'.format(self.chainID1,self.resseq1,self.chainID2,self.resseq2)
                pass
            else:
                return '### patch unknown for '+str(self)+'\n'


