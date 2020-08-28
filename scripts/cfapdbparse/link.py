from segment import _seg_class_
from residue import get_residue,get_atom
class Link:
    def __init__(self,pdbrecord=[]):
        self.pdbrecord=pdbrecord
# 1 -  6         Record name    "LINK  "
        if len(pdbrecord)>0:
            self.record_name=pdbrecord[0:6].strip()
#13 - 16         Atom           name1           Atom name.
            self.name1=pdbrecord[12:16].strip()
#17              Character      altLoc1         Alternate location indicator.
            self.altloc1=pdbrecord[16:17]
#18 - 20         Residue name   resName1        Residue  name. MODIFIED for 4-char resnames
            self.resname1=pdbrecord[17:21].strip()
#    22              Character      chainID1        Chain identifier.
            self.chainID1=pdbrecord[21:22]
#23 - 26         Integer        resSeq1         Residue sequence number.
            self.resseqnum1=int(pdbrecord[22:26])
#27              AChar          iCode1          Insertion code.
            self.icode1=pdbrecord[26:27]
#43 - 46         Atom           name2           Atom name.
            self.name2=pdbrecord[42:46].strip()
#47              Character      altLoc2         Alternate location indicator.
            self.altloc2=pdbrecord[46:47]
#48 - 50         Residue name   resName2        Residue name.  MODIFIED for 4-char resnames
            self.resname2=pdbrecord[47:51].strip()
#52              Character      chainID2        Chain identifier.
            self.chainID2=pdbrecord[51:52]
#53 - 56         Integer        resSeq2         Residue sequence number.
            self.resseqnum2=int(pdbrecord[52:56])
#57              AChar          iCode2          Insertion code.
            self.icode2=pdbrecord[56:57]
#60 - 65         SymOP          sym1            Symmetry operator atom 1.
            self.sym1=pdbrecord[59:65].strip()
#67 - 72         SymOP          sym2            Symmetry operator atom 2.
            self.sym2=pdbrecord[66:72].strip()
#74 – 78         Real(5.2)      Length          Link distance
            self.link_distance=float(pdbrecord[73:78])
            self.segname1=self.chainID1
            self.segname2=self.chainID2
            self.residue1=''
            self.residue2=''
            self.atom1=''
            self.atom2=''
            self.empty=False
        else:
            self.empty=True
    def pdb_line(self):
        pdbline='{:6s}'.format(self.record_name)+6*' '+'{:>4s}'.format(self.name1+' ' if len(self.name1)<3 else self.name1)+'{:1s}'.format(self.altloc1)+'{:3s}'.format(self.resname1)+' '+'{:1s}'.format(self.chainID1)+'{:4d}'.format(self.resseqnum1)+'{:1s}'.format(self.icode1)+16*' '+'{:4>s}'.format(self.name2+' ' if len(self.name2)<3 else self.name2)+'{:1s}'.format(self.altloc2)+'{:3s}'.format(self.resname2)+' '+'{:1s}'.format(self.chainID2)+'{:4d}'.format(self.resseqnum2)+'{:1s}'.format(self.icode2)+2*' '+'{:>6s}'.format(self.sym1)+' '+'{:>6s}'.format(self.sym2)+'{:6.2f}'.format(self.link_distance)
        return pdbline
    def updateSegnames(self,R):
        self.residue1=get_residue(R,self.chainID1[0],self.resseqnum1)
        self.residue2=get_residue(R,self.chainID2[0],self.resseqnum2)
        self.atom1=get_atom(R,self.chainID1[0],self.resseqnum1,self.name1)
        self.atom2=get_atom(R,self.chainID2[0],self.resseqnum2,self.name2)
        self.segname1=self.atom1.segname
        self.segname2=self.atom2.segname
    def isInLink(self,chain,resid,pos=''):
        if pos=='':
            if (self.chainID1==chain and self.resseqnum1==resid) or (self.chainID2==chain and self.resseqnum2==resid):
               return True
            else:
               return False
        elif pos==1:
            if self.chainID1==chain and self.resseqnum1==resid:
               return True
            else:
               return False
        elif pos==2:
            if self.chainID2==chain and self.resseqnum2==resid:
               return True
            else:
               return False
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
        return retstr.format(self.record_name,self.name1,self.altloc1,self.resname1,self.chainID1,self.resseqnum1,self.icode1,self.name2,self.altloc2,self.resname2,self.chainID2,self.resseqnum2,self.icode2,self.sym1,self.sym2,self.link_distance)
    def psfgen_patchline(self):
        if self.resname1=='ASN' and _seg_class_[self.resname2]=='GLYCAN':
            return 'patch NGLB {}:{} {}:{}\n'.format(self.segname1,self.resseqnum1,self.segname2,self.resseqnum2)
        else:
            # for a glycan-glycan patch, the C1 atom is always on the ji-residue
            if self.name2=='C1' and _seg_class_[self.resname1]=='GLYCAN':
                cmdi='[axeq {} 0 {} {} {}]'.format(self.resseqnum2,self.chainID2,self.name2,self.resseqnum1)
                cmdj='[axeq {} 0 {} {} {}]'.format(self.resseqnum1,self.chainID1,self.name1,-1)
                return 'patch 1{:1s}{}{} {}:{} {}:{}\n'.format(self.name1[1], cmdi,cmdj,self.segname1,self.resseqnum1,self.segname2,self.resseqnum2)
            elif self.name1=='C1' and _seg_class_[self.resname2]=='GLYCAN':
                cmdj='[axeq {} 0 {} {} {}]'.format(self.resseqnum2,self.chainID2,self.name2,self.resseqnum1)
                cmdi='[axeq {} 0 {} {} {}]'.format(self.resseqnum1,self.chainID1,self.name1,-1)           
                return 'patch 1{:1s}{}{} {}:{} {}:{}\n'.format(self.name2[1], cmdi,cmdj,self.segname2,self.resseqnum2,self.segname1,self.resseqnum1)
            elif self.name1=='O6' and self.name2=='C2':
                return 'patch SA26E {}:{} {}:{}\n'.format(self.segname1,self.resseqnum1,self.segname2,self.resseqnum2)
                pass
            else:
                return '### patch unknown for '+str(self)+'\n'


