import operator
from segment import Segment, _seg_typedict_byresname_, _segname_second_character_
class Chain:
    def __init__(self,r=None,parent_molecule=None,source_chainID='',cleavage_chainID=''):
        ''' source_chainID is this chain's chain ID in an input PDB file 
            cleavage_chainID is this chain's local chainID prior to 
            any cleavage. '''
        self.residues=[]
        self.Segments=[]
        self.subCounter={}
        self.subCounter['GLYCAN']=0
        self.subCounter['LIGAND']=0
        self.subCounter['PROTEIN']=0
        self.subCounter['ION']=0
        self.subCounter['WATER']=0
        if r!=None:
            self.chainID=r.chainID
            self.residues=[r]
            if source_chainID!=None:
                self.source_chainID=source_chainID
            else:
                self.source_chainID=r.chainID
        self.parent_molecule=parent_molecule
        if cleavage_chainID!=None:
            self.cleavage_chainID=cleavage_chainID
        else:
            self.cleavage_chainID=self.chainID

    def get_molid(self):
        return self.parent_molecule.molid

    def add_residue(self,r):
        if len(self.residues)==0:
            self.residues=[r]
            self.chainID=r.chainID
        elif self.chainID==r.chainID:
            self.residues.append(r)
    def sort_residues(self):
        ''' sort list of residues by resseqnum '''
        sorted_residues=sorted(self.residues,key=operator.attrgetter('resseqnum','insertion'))
        self.residues=sorted_residues
        mn=99999
        mx=-99999
        mnr=''
        mxr=''
        self.Nterm=''
        self.Cterm=''
        for r in self.residues:
            if _seg_typedict_byresname_[r.name]=='PROTEIN':
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
            if _seg_typedict_byresname_[r.name]=='PROTEIN' and len(r.down)>0:
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
            if _seg_typedict_byresname_[r.name]=='PROTEIN':
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
    def MakeSegments(self,Links):
        ''' scans residues in a chain to divvy them up into segments.
            A chain can comprise one or more distinct segments.  It is assumed
            however that there is only one PROTEIN segment in a chain, and its
            segment name is its chain name.  Other segments belonging to the chain 
            could be IONs, WATERs, LIGANDs, or GLYCANs, and their segment names
            are derived from the parent chain names. '''
        Segments=[]
        for r in self.residues:
            thissegtype=_seg_typedict_byresname_[r.name]
            if Segments==[]:
                Segments.append(Segment(r,subcounter=self.nextSubCounter(thissegtype),parent_chain=self))
            else:
                ''' find the segment to add this residue to '''
                for s in Segments:
                    found=False
                    ''' '''
                    if s.segtype==thissegtype:
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
                    # initiate a new segment
                    Segments.append(Segment(r,subcounter=self.nextSubCounter(thissegtype),parent_chain=self))
            #print('#### MakeSegments: residue {} placed in segment {}'.format(str(r),str(s)))
        # apportion mutations in this chain to their correct segments
        for m in self.Mutations:
            found=False
            for s in Segments:
                if s.segtype=='PROTEIN':
                    for r in s.residues:
                        if r.resseqnum==m.resseqnum and r.insertion==m.insertion:
                            s.mutations.append(m)
                            found=True
            if not found:
                print('#### Warning: no protein segment in chain {} found for mutation {}'.format(self.chainID,str(m)))
        # apportion deletions to their correct segments
        for d in self.Deletions:
            found=False
            for s in Segments:
                if s.segtype=='PROTEIN':
                    for r in s.residues:
                        if r.resseqnum==d.resseqnum and r.insertion==m.insertion:
                            s.deletions.append(d)
                            found=True
            if not found:
                print('#### Warning: no protein segment in chain {} found for deletion {}'.format(self.chainID,str(d)))
        # apportion grafts to their correct segments
        for g in self.Grafts:
            found=False
            g.target_segment=None
            for s in self.Segments:
                for r in s.residues:
                    if r.resseqnum==g.target_res and r.insertion==g.target_ins:
                        g.target_segment=s
                        s.apply_graft(g)
                        found=True
            if not found:
                print('#### Warning: no target residue {} found in segments of chain {} for graft {}'.format(g.target_res,g.target_chain,str(g)))
        # create a new segment for each attachment
#        for a in Attachments:
 #           self.ImportAttachmentSegment(a)
        return Segments
        
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
        retstr=''
        if segtype!='PROTEIN':
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

if __name__=='__main__':
    # testing
    class Mock:
        def __init__(self,rsn=0,ic='',rn='',chainID='A'):
            self.resseqnum=rsn
            self.name=rn
            self.insertion=ic
            self.chainID=chainID
        def __str__(self):
            return '{}_{}{}{}'.format(self.chainID,self.name,self.resseqnum,self.insertion)
    someresids=[101,103,102,105,101,105,106,101]
    someics   =['A','' ,'' ,'' ,'' ,'A','' ,'B']
    C=Chain()
    for r,i in zip(someresids,someics):
        C.add_residue(Mock(rsn=r,ic=i,rn='ALA'))
    print('Unsorted:')
    for r in C.residues:
        print(str(r))
    C.sort_residues()
    print('Sorted:')
    for r in C.residues:
        print(str(r))

    