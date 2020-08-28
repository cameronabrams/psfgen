import operator
from segment import Segment,_seg_class_

class Chain:
    def __init__(self,r):
        self.chainID=r.chainID
        self.residues=[r]
        self.source_chainID=''
        self.Segments=[]
        self.subCounter={}
        self.subCounter['GLYCAN']=0
    def add_residue(self,r):
        if self.chainID==r.chainID:
            self.residues.append(r)
    def sort_residues(self):
        ''' sort list of residues by resseqnum '''
        sorted_residues=sorted(self.residues,key=operator.attrgetter('resseqnum'))
        self.residues=sorted_residues
        mn=99999
        mx=-99999
        mnr=''
        mxr=''
        self.Nterm=''
        self.Cterm=''
        for r in self.residues:
            if _seg_class_[r.name]=='PROTEIN':
                if r.resseqnum<mn:
                    mn=r.resseqnum
                    mnr=r
                if r.resseqnum>mx:
                    mx=r.resseqnum
                    mxr=r
#        if mnr=='' or mxr=='':
#            print('Note, no protein residues found chain {}'.format(self.chainID))
        self.Nterm=mnr
        self.Cterm=mxr
    def group_residues(self):
        ''' group residues by connectivity '''
        owners=[]
        for i,r in enumerate(self.residues):
            if _seg_class_[r.name]=='PROTEIN' and len(r.down)>0:
                owners.append([r,i])
        for oi in owners:
            o,i=oi
            #print(o)
            #print('{}{}-{}{}->'.format(i,o.chainID,o.name,o.resseqnum))
            for dd in o.get_down_group():
                if dd.chainID==o.chainID:
                #print('    {}-{}{}'.format(dd.chainID,dd.name,dd.resseqnum))
                    self.residues.remove(dd)
                    self.residues.insert(i+1,dd)
                    i+=1
    def __str__(self):
        protein_resid=[]
        other_resid=[]
        for r in self.residues:
            if _seg_class_[r.name]=='PROTEIN':
                protein_resid.append(r.resseqnum)
            else:
                other_resid.append(r.resseqnum)
        rstr='CHAIN {} {:d} - {:d}'.format(self.chainID,min(protein_resid),max(protein_resid))
        if len(other_resid)>0:
            rstr+=' + '
            for o in other_resid:
                rstr+='{:d}{}'.format(o,', ' if o!=other_resid[-1] else '')
        return rstr
    def Cleave(self,Cleave):
        C_index=-1
        self.group_residues()
        for i,r in enumerate(self.residues):
            if r.resseqnum==Cleave.parent_Cterm_resseqnum:
                C_index=i
        dres=self.residues[C_index+1:]
        self.residues=self.residues[0:C_index+1]
        for r in dres:
            r.set_chainID(Cleave.daughter_chainID)
        Daughter=Chain(dres[0])
        Daughter.residues.extend(dres[1:])
        self.sort_residues()
        Daughter.sort_residues()
        Daughter.source_chainID=self.chainID
        return Daughter
    def MakeSegments(self,Links,Mutations=[],Grafts=[]):
        self.Segments=[]
        for r in self.residues:
            if self.Segments==[]:
                if _seg_class_[r.name]=='GLYCAN':
                   s=Segment(r,self.nextSubCounter('GLYCAN'),chain=self)
                else:
                   s=Segment(r,chain=self)
                self.Segments.append(s)
            else:
                for s in self.Segments:
                    #fp.write('looking {} {}'.format(s.segtype,_seg_class_[r.name]))
                    found=False
                    # if this is not an ion and not a water, then we will add this residue to the 
                    # segment that contains a residue that is covalently linked to it via a LINK
                    # ions and waters are each added to a single segment per chain
                    # glycans are added to separate segments
                    if s.segtype==_seg_class_[r.name]:
                        if s.segtype=='GLYCAN':
                            # only add this residue to this segment if it is linked to an
                            # existing residue in that segment 
                            if self.isConnected(r,s,Links):
                                s.add_residue(r)
                                found=True
                                break
                        else:
                            #fp.write('adding')
                            s.add_residue(r)
                            found=True
                            break
                if not found:
                    # if this is an ion, segname is AI, where A is chainID, 'I' means ion
                    # if this is a water, segname is AW
                    # if this is a glycan, segname is AGnn, where nn is the next avail. glycan number
                    if _seg_class_[r.name]=='GLYCAN':
                       s=Segment(r,self.nextSubCounter('GLYCAN'),chain=self)
                    else:
                       s=Segment(r,chain=self)
                    self.Segments.append(s)
        if len(Mutations)>0:
            for m in Mutations:
                found=False
                if m.chainID==self.chainID:
                    for s in self.Segments:
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
        if len(Grafts)>0:
            for g in Grafts:
                g.target_segment=''
                found=False
                if g.target_chain==self.chainID:
                    for s in self.Segments:
                        for r in s.residues:
                            if r.resseqnum==g.target_res:
                                found=True
                                g.target_segment=s
                                s.apply_graft(g)
                                break
                            else:
                                pass
                else:
                    pass
        else:
            pass

    def nextSubCounter(self,segtype):
        retstr='{:02d}'.format(self.subCounter[segtype])
        self.subCounter[segtype]+=1
        return retstr
    def isConnected(self,r,s,l):
        con=False
        for jr in s.residues:
            if r.chainID==jr.chainID:
                for li in l:
                    if li.chainID1==r.chainID:
                        if (r.resseqnum==li.resseqnum1 and jr.resseqnum==li.resseqnum2) or (r.resseqnum==li.resseqnum2 and jr.resseqnum==li.resseqnum1):
                           # print(r,'connected to',jr,'via',li)
                            con=True
                            break
            else:
                pass
        return con
    def has_resseqnum(self,resseqnum):
        for r in self.residues:
            if r.resseqnum==resseqnum:
                return True
        return False
