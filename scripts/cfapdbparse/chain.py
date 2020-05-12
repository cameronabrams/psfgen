import operator
from segment import Segment,_seg_class_

class Chain:
    def __init__(self,r):
        self.chainID=r.chainID
        self.residues=[r]
        self.source_chainID=r.chainID
    def add_residue(self,r):
        if self.chainID==r.chainID:
            self.residues.append(r)
    def sort_residues(self):
        sorted_residues=sorted(self.residues,key=operator.attrgetter('resseqnum'))
        self.residues=sorted_residues
    def __str__(self):
        resid=[]
        for r in self.residues:
            resid.append(r.resseqnum)
        return 'CHAIN {} {:d} - {:d}'.format(self.chainID,min(resid),max(resid))
    def Cleave(self,Cleave):
        C_index=-1
        for i,r in enumerate(self.residues):
            if r.resseqnum==Cleave.parent_Cterm_resseqnum:
                C_index=i
        dres=self.residues[C_index+1:]
        self.residues=self.residues[0:C_index+1]
        for r in dres:
            r.set_chainID(Cleave.daughter_chainID)
        Daughter=Chain(dres[0])
        Daughter.residues.extend(dres[1:])
        return Daughter
    def Segments(self,Mutations=0):
        S=[]
        for r in self.residues:
            if S==[]:
                S.append(Segment(r))
            else:
                for s in S:
                    #fp.write('looking {} {}'.format(s.segtype,_seg_class_[r.name]))
                    found=False
                    if s.segtype==_seg_class_[r.name]:
                        #fp.write('adding')
                        s.add_residue(r)
                        found=True
                        break
                if not found:
                    #fp.write('creating')
                    S.append(Segment(r))
        if Mutations != 0:
            for m in Mutations:
                found=False
                if m.chainID==self.chainID:
                    for s in S:
                        if s.segtype=='PROTEIN':
                            for r in s.residues:
                                if r.resseqnum==m.resseqnum:
                                    found=True
                                    s.mutations.append(m)
                                    break
                                else:
                                    pass
                        else:
                            pass
                else:
                    pass
        else:
            pass

        return S


