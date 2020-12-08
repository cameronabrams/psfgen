import operator
from segment import Segment,_seg_class_,_segname_second_character_
class Chain:
    def __init__(self,r,parent_molecule=''):
        self.chainID=r.chainID
        self.residues=[r]
        self.source_chainID=r.chainID # has a value if C-term-product of cleavage
        self.Segments=[]
        self.subCounter={}
        self.subCounter['GLYCAN']=0
        self.subCounter['LIGAND']=0
        self.parent_molecule=parent_molecule

    def get_molid(self):
        return self.parent_molecule.molid

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
        Daughter.parent_molecule=self.parent_molecule
        return Daughter
    def MakeSegments(self,Links,Mutations=[],Grafts=[],Attachments=[]):
        ''' scans residues in a chain to divvy them up into segments '''
        self.Segments=[]
        for r in self.residues:
            if self.Segments==[]:
                ''' create the first segment '''
                if _seg_class_[r.name]=='GLYCAN':
                   ''' if this is a glycan residue, use the segnaming convention '''
                   s=Segment(r,subcounter=self.nextSubCounter('GLYCAN'),parent_chain=self)
                elif _seg_class_[r.name]=='LIGAND':
                   ''' if this is a ligand residue, use the segnaming convention '''
                   s=Segment(r,subcounter=self.nextSubCounter('LIGAND'),parent_chain=self)
                else:
                   s=Segment(r,parent_chain=self)
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
                        elif s.segtype!='LIGAND':
                            #fp.write('adding')
                            s.add_residue(r)
                            found=True
                            break
                if not found:
                    # if this is an ion, segname is AI, where A is chainID, 'I' means ion
                    # if this is a water, segname is AW
                    # if this is a glycan, segname is AGnn, where nn is the next avail. glycan number
                    if _seg_class_[r.name]=='GLYCAN':
                       s=Segment(r,subcounter=self.nextSubCounter('GLYCAN'),parent_chain=self)
                    elif _seg_class_[r.name]=='LIGAND':
                       s=Segment(r,subcounter=self.nextSubCounter('LIGAND'),parent_chain=self)
                    else:
                       s=Segment(r,parent_chain=self)
                    self.Segments.append(s)
            #print('#### MakeSegments: residue {} placed in segment {}'.format(str(r),str(s)))
        # apportion mutations to their correct segments
        for m in [_ for _ in Mutations if _.chainID==self.chainID]:
            found=False
            for s in self.Segments:
                if s.segtype=='PROTEIN':
                    for r in s.residues:
                        if r.resseqnum==m.resseqnum:
                            s.mutations.append(m)
                            found=True
            if not found:
                print('#### Warning: no protein segment in chain {} found for mutation {}'.format(self.chainID,str(m)))
        # apportion grafts to their correct segments
        for g in [_ for _ in Grafts if _.target_chain==self.chainID]:
            found=False
            g.target_segment=''
            for s in self.Segments:
                for r in s.residues:
                    if r.resseqnum==g.target_res:
                        g.target_segment=s
                        s.apply_graft(g)
                        found=True
            if not found:
                print('#### Warning: no target residue {} found in segments of chain {} for graft {}'.format(g.target_res,c.chainID,str(g)))
        # create a new segment for each attachment
        for a in Attachments:
            self.ImportAttachmentSegment(a)

    def ImportAttachmentSegment(self,a):
        ''' an attachment segment in a chain in an auxiliary molecule is imported
            into chain self '''
        newseg=a.source_segment
        ''' detch from its own molecule '''
        newseg.parent_chain=self
        if newseg.segtype=='GLYCAN':
            newseg.segname=newseg.get_chainID()+_segname_second_character_['GLYCAN']+self.nextSubCounter('GLYCAN')
        else:
            newseg.segname=newseg.get_chainID()
        for r in newseg.residues:
            r.segname=newseg.segname
            r.chainID=newseg.get_chainID()
            for a in r.atoms:
                a.segname=newseg.segname
                a.chainID=newseg.get_chainID()
        self.attach=a
        self.Segments.append(newseg)

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
    def relabel(self,newchainID):
        self.chainID=newchainID
        for r in self.residues:
            r.chainID=newchainID
            for a in r.atoms:
                a.chainID=newchainID

