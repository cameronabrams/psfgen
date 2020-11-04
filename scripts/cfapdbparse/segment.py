from residue import  ResnameCharmify,_PDBResName123_, _pdb_glycans_, _pdb_ions_, _ResNameDict_PDB_to_CHARMM_, _ResNameDict_CHARMM_to_PDB_
_seg_class_={'HOH':'WATER'}
_seg_class_.update({k:'ION' for k in _pdb_ions_})
_seg_class_.update({k:'GLYCAN' for k in _pdb_glycans_})
_seg_class_.update({_ResNameDict_PDB_to_CHARMM_[k]:'GLYCAN' for k in _pdb_glycans_})
_seg_class_.update({k:'PROTEIN' for k in _PDBResName123_.values()})
_seg_class_['HIS']='PROTEIN'
_seg_class_['HSE']='PROTEIN'
_seg_class_['HSD']='PROTEIN'
_segname_second_character_={'PROTEIN':'','ION':'I','WATER':'W','GLYCAN':'G','OTHER':'O'}
import sel

class SubsegmentBounds:
    def __init__(self,l=-1,r=-1,typ='NONE',d=''):
        self.l=l
        self.r=r
        self.typ=typ
        self.d=d

class Fragment:
    ''' a set of contiguous residues with no gaps that can be loaded into a psfgen
        segment stanza by invoking a pdb file 
    '''
    def __init__(self,source_chainID,replica_chainID,resseqnum1,resseqnum2):
        self.source_chainID=source_chainID
        self.replica_chainID=replica_chainID
        self.resseqnum1=resseqnum1
        self.resseqnum2=resseqnum2
    def pdb_str(self):
        return '{}_{}_to_{}.pdb'.format(self.replica_chainID,self.resseqnum1,self.resseqnum2)
    def __str__(self):
        return 'FRAGMENT: {} {} to {}'.format(self.replica_chainID,self.resseqnum1,self.resseqnum2)

class Loop:
    ''' a set of contiguous residues from REMARK 465 pdb entries; i.e., they
        are missing from the base pdb file but present in the construct.  They
        must be included in a psfgen segment stanza via the 'residue' invocation.
    '''
    def __init__(self,source_chainID,replica_chainID,resseqnum0,r):
        self.source_chainID=source_chainID
        self.replica_chainID=replica_chainID
        self.resseqnum0=resseqnum0
        self.residues=[r]
        self.term='UNSET'
    def add_residue(self,r):
        self.residues.append(r)
    def __str__(self):
        return 'LOOP: {} ({})-[{} to {}]'.format(self.replica_chainID,self.resseqnum0,self.residues[0].resseqnum,self.residues[-1].resseqnum)
    def caco_str(self):
        return 'coord {} {} N [cacoIn_nOut {} {} 0]\n'.format(self.replica_chainID,self.residues[0].resseqnum,self.resseqnum0,self.replica_chainID)

class Segment:
    """A class for holding all information necessary to generate a segment stanza in psfgen.

       Class instance attributes:
       * segname: name of segment
       * parent_chain: chain to which segment belongs
       * segtype: 'PROTEIN' or 'GLYCAN', determined by _seg_class_ global
       * residues: list of residues
       * mutations: list of mutations
       * graft: graft designation; if set, will contain instructions on how to build segment from a graft
       
       """
    def __init__(self,r,subcounter='',parent_chain=''):
        """Initializes a segment instance by passing in first residue of segment"""
        self.segname=r.chainID+_segname_second_character_[_seg_class_[r.name]]+subcounter
        self.parent_chain=parent_chain
        self.segtype=_seg_class_[r.name]
        self.residues=[r]
        self.mutations=[]
        self.graft=''
        self.rootres=''
        self.attach=''
        self.pdbfiles=[]
        if _seg_class_[r.name]=='GLYCAN':
            #print(self.segname,r.name,r.resseqnum,r.up)
            self.rootres=r.up[0]
        r.segname=self.segname
        for a in r.atoms:
            a.segname=self.segname
    def get_chainID(self):
        return self.parent_chain.source_chainID
    def get_molecule(self):
#        print(self.parent_chain.parent_molecule)
        return self.parent_chain.parent_molecule
    def __str__(self):
        retbase='({}){}[{}] {} - {}'.format(self.parent_chain.chainID,self.segname,self.segtype,self.residues[0].resseqnum,self.residues[-1].resseqnum)
        if self.rootres!='':
            r=self.rootres
            retbase='('+str(r)+')->'+retbase
        if self.graft!='':
            retbase+='::'+self.graft.source_pdb+':'+str(self.graft.source_segment)
        return retbase
    def add_residue(self,r):
        self.residues.append(r)
        r.segname=self.segname
        for a in r.atoms:
            a.segname=self.segname
    def apply_graft(self,g):
        self.graft=g
    def subsegments(self,source_chainID,replica_chainID):
        self.subsegbounds=[]
        Loops=[]
        curr=SubsegmentBounds() # bounds as indices in self.residues[]
        for i,r in enumerate(self.residues):
             if len(r.atoms)>0:
                 if curr.typ=='NONE':
                     curr.typ='FRAGMENT'
                     curr.l=i
                     curr.r=i
                 elif curr.typ=='LOOP':
                     self.subsegbounds.append(curr)
                     curr=SubsegmentBounds(i,i,'FRAGMENT','NULL')
                 elif curr.typ=='FRAGMENT':
                     curr.r=i
             else:
                 if curr.typ=='NONE':
                     curr.typ='LOOP'
                     curr.l=i
                     curr.r=i
                 elif curr.typ=='FRAGMENT':
                     self.subsegbounds.append(curr)
                     curr=SubsegmentBounds(i,i,'LOOP','NULL')
                 elif curr.typ=='LOOP':
                     curr.r=i
        self.subsegbounds.append(curr)
        lst=-1
        for j,b in enumerate(self.subsegbounds):
             if b.typ=='FRAGMENT':
                 b.d=Fragment(source_chainID,replica_chainID,self.residues[b.l].resseqnum,self.residues[b.r].resseqnum)
                 lst=b.r
             elif b.typ=='LOOP':
                 L=Loop(source_chainID,replica_chainID,self.residues[lst].resseqnum,self.residues[b.l])
                 for i in range(b.l+1,b.r+1):
                     L.add_residue(self.residues[i])
                 Loops.append(L)
                 b.d=L
                 if j==0 or j==len(self.subsegbounds)-1:
                     b.d.term=False # not a terminated loop (this is an end!)
                 else:
                     b.d.term=True
                 #print('{} {} {} {} {}'.format(j,b.d.chainID,b.d.residues[0].resseqnum,b.d.residues[-1].resseqnum,'terminated' if b.d.term else 'not terminated'))
        return Loops

    def write_psfgen_stanza(self,includeTerminalLoops=False,tmat=None):
        stanzastr=''
        my_chainID=self.get_chainID()
        rep_chainID=tmat.get_replica_chainID(my_chainID)
        rep_segname=self.segname.replace(my_chainID,rep_chainID,1)
        #print('#### writing stanza for chain {} (source {}) segname {}'.format(rep_chainID,my_chainID,rep_segname))
        if tmat==None:
            print('ERROR: write_psfgen_stanza needs a tmat!')
            exit()
        if self.segtype=='PROTEIN':
            if self.graft!='' or self.attach!='':
                print('ERROR: Protein grafts and attachment segments are not yet implemented')
                return 'ERROR',[]
            Loops=self.subsegments(my_chainID,rep_chainID)
            ''' PART 1:  Process selections and write PDB files for fragments '''
            for ss in self.subsegbounds:
                if ss.typ=='FRAGMENT':
                    r=ss.d
                    stanzastr+='set mysel [atomselect ${} "chain {} and resid {} to {}"]\n'.format(self.get_molecule().molid_varname,self.parent_chain.source_chainID,r.resseqnum1,r.resseqnum2)
                    if not tmat.isidentity():
                         stanzastr+=sel.backup('mysel')
                         stanzastr+='$mysel move {}\n'.format(tmat.OneLiner())
                    ''' tmat transformation '''
                    stanzastr+='$mysel writepdb {}\n'.format(r.pdb_str())
                    if not tmat.isidentity():
                         stanzastr+=sel.restore('mysel')
                    self.pdbfiles.append(r.pdb_str())
            ''' PART 2:  Build segment stanza '''
            stanzastr+='segment {} {{\n'.format(rep_segname)
            for i,ss in enumerate(self.subsegbounds):
                if ss.typ=='FRAGMENT':
                    f=ss.d
                    stanzastr+='   pdb {}\n'.format(f.pdb_str())
                elif ss.typ=='LOOP':
                    if (i==0 or i==(len(self.subsegbounds)-1)) and not includeTerminalLoops:
                       ''' shunt this if this is a terminal loop and includeTerminalLoops is False '''
                       pass
                    else:
                       ''' this is either NOT a terminal loop, or if it is, includeTerminalLoops is True '''
                       l=ss.d
                       for rr in l.residues:
                           nm=ResnameCharmify(rr.name)
                           stanzastr+='   residue {}{} {} {}\n'.format(rr.resseqnum,rr.insertion,nm,tmat.get_replica_chainID(rr.chainID))
            ''' PART 2.1:  Include mutations '''
            #print('### {} mutations'.format(len(self.mutations)))
            for m in self.mutations:
                stanzastr+=m.psfgen_segment_str()
            stanzastr+='}\n'
            ''' PART 3:  Issue coordinate-setting commands '''
            for i,ss in enumerate(self.subsegbounds):
                if ss.typ=='FRAGMENT':
                    stanzastr+='coordpdb {} {}\n'.format(ss.d.pdb_str(),rep_segname)
                elif ss.typ=='LOOP':
                    if (i==0 or i==(len(self.subsegbounds)-1)) and not includeTerminalLoops:
                        pass
                    else:
                        stanzastr+=ss.d.caco_str()
            return stanzastr,Loops
        elif self.segtype=='GLYCAN':
            stanzastr=''
            if self.graft!='':
                ''' this is a segment with an associated graft '''
                g=self.graft
                g.ingraft_segname=self.segname # graft inherits segname
                g.ingraft_chainID=self.parent_chain.chainID
                stanzastr+=g.transform(self.parent_chain.parent_molecule)
                self.pdbfiles.append(g.transformed_pdb)
                ''' g.transformed_pdb is now the file with inputs! '''
                stanzastr+=r'segment {} {{'.format(rep_segname)+'\n'
                stanzastr+='     pdb {}\n'.format(g.transformed_pdb)
                stanzastr+='}\n'
                stanzastr+='coordpdb {} {}\n'.format(g.transformed_pdb,rep_segname)
            elif self.attach!='':
                ''' this is a segment from another molecule to be attached '''
                pass
            else:
                ''' this is an existing glycan segment '''
                f=Fragment(my_chainID,rep_chainID,self.residues[0].resseqnum,self.residues[-1].resseqnum)
                pdb=f.pdb_str()
                self.pdbfiles.append(pdb)
                #print('#### segment {}'.format(self.segname))
                stanzastr+='set mysel [atomselect ${} "chain {} and resid {} to {}"]\n'.format(self.get_molecule().molid_varname,self.parent_chain.source_chainID,self.residues[0].resseqnum,self.residues[-1].resseqnum)
                if not tmat.isidentity():
                     stanzastr+=sel.backup('mysel')
                     stanzastr+='$mysel move {}\n'.format(tmat.OneLiner())
                stanzastr+=sel.charmm_namify('mysel')
                stanzastr+='$mysel writepdb {}\n'.format(pdb)
                if not tmat.isidentity():
                     stanzastr+=sel.restore('mysel')
                stanzastr+='segment {} {{\n'.format(rep_segname)
                stanzastr+='    pdb {}\n'.format(pdb)
                stanzastr+='}\n'
                stanzastr+='coordpdb {} {}\n'.format(pdb,rep_segname)
            return stanzastr,[]
        elif self.segtype in ['WATER','ION', 'OTHER']:
            f=Fragment(my_chainID,tmat.get_replica_chainID(my_chainID),self.residues[0].resseqnum,self.residues[-1].resseqnum)
            pdb=f.pdb_str()
            self.pdbfiles.append(pdb)
            stanzastr+='set mysel [atomselect ${} "chain {} and resid {} to {}"]\n'.format(self.get_molid(),my_chainID,self.residues[0].resseqnum,self.residues[-1].resseqnum)
            stanzastr+=sel.charmm_namify('mysel')
            if not tmat.isidentity():
                 stanzastr+=sel.backup('mysel')
                 stanzastr+='$mysel move {}\n'.format(tmat.OneLiner())
            stanzastr+='$mysel writepdb {}\n'.format(pdb)
            if not tmat.isidentity():
                 stanzastr+=sel.restor('mysel')
            stanzastr+='segment {} {{\n'.format(rep_segname)
            stanzastr+='    pdb {}\n'.format(thispdb)
            stanzastr+='}\n'
            stanzastr+='coordpdb {} {}\n'.format(thispdb,rep_segname)
            return stanzastr,[] 
        else:
            print('ERROR: Unrecognized segment type {}'.format(self.segtype))
            return 'ERROR',[]

