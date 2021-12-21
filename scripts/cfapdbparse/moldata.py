from os import kill
from atom import Atom, _PDBAtomNameDict_
from ssbond import SSBond
from missing import Missing
from link import Link
from chain import Chain
from seqadv import Seqadv
from mutation import Mutation
from residue import Residue, _PDBResName123_, _pdb_glycans_, _pdb_ions_, _ResNameDict_PDB_to_CHARMM_, _ResNameDict_CHARMM_to_PDB_, get_residue, get_atom
from segment import Segment, _seg_typedict_byresname_

def _recursiveDelDown(R,L,r):
    for rd in r.down:
        L.remove(r.downlink[r.down.index(rd)])
        _recursiveDelDown(R,L,rd)
    print('Deleting residue',r.printshort())
    R.remove(r)

class MolData:
    ''' data in a pdb/cif file that is inheritable by biological assemblies, e.g.,
        atoms, residues, chains, links, ssbonds '''
    PDBClonables=['ATOM','HETATM','LINK','SSBOND','SEQADV']
    def __init__(self,userMods={},parent_molecule=''):
        self.userMods=userMods
        self.Atoms=[]
        self.Residues=[]
        self.Chains={} # keyed by chain id 'A', 'B', 'C', etc.
        self.Seqadv=[] # Sequence variations from databases explicitly listed in PDB file
        self.Links=[] if not 'userLinks' in userMods else userMods['userLinks']
        self.SSBonds=[] if not 'userSSBonds' in userMods else userMods['userSSBonds']
        self.MissingRes=[] if not 'userMissing' in userMods else userMods['userMissing']
        self.Mutations=[] if not 'userMutations' in userMods else userMods['userMutations']
        self.Deletions=[] if not 'userDeletions' in userMods else userMods['userDeletions']
        self.Grafts=[] if not 'userGrafts' in userMods else userMods['userGrafts']
        self.Cleavages=[] if not 'userCleavages' in userMods else userMods['userCleavages']
        self.userXSSBonds=[] if not 'userXSSBonds' in userMods else userMods['userXSSBonds']
        self.parent_molecule=parent_molecule
        self.adjustedChainIDs=[]
    def Parse(self,line):
        key=line[:6].strip()
        if key=='ATOM' or key=='HETATM':
            self.Atoms.append(Atom(line))
        elif key=='LINK':
            self.Links.append(Link(line))
        elif key=='SSBOND':
            self.SSBonds.append(SSBond(line))
        elif key=='SEQADV':
            self.Seqadv.append(Seqadv(line))
        elif key=='REMARK':
            if line[7:10]=='465':
                test_int=line[20:26].strip()
                if test_int.isdigit() or (len(test_int)>0 and test_int[0]=='-'):
                    self.MissingRes.append(Missing(line))
    def MakeConstructs(self):
        self._MakeResidues()
        fixConflicts=False if not 'fixConflicts' in self.userMods else self.userMods['fixConflicts']
        fixEngineeredMutations=False if not 'fixEngineeredMutations' in self.userMods else self.userMods['fixEngineeredMutations']
        for sa in self.Seqadv:
            if fixConflicts==True and sa.conflict=='CONFLICT':
#                print('adding as mutation: {}'.format(sa.printshort()))
                self.Mutations.append(Mutation(seqadv=sa))
            elif fixEngineeredMutations==True and sa.conflict=='ENGINEERED MUTATION':
                self.Mutations.append(Mutation(seqadv=sa)) 
#                print('adding as mutation: {}'.format(sa.printshort()))
        self._MakeChains()
        self._MakeLinks()
        self._checkMutationsAgainstLinks()
        self._checkDeletionsAgainstLinks()
    def ShowSeqadv(self,brief=False,indent=''):
        if len(self.Seqadv)>0:
            if not brief:
                for sa in self.Seqadv:
                    print('### SEQADV:',sa)
            else:
                counts={}
                for sa in self.Seqadv:
                    if sa.conflict not in counts:
                        counts[sa.conflict]=1
                    else:
                        counts[sa.conflict]+=1
                print(f'{indent}{sum(counts.values())} SEQADV records:')
                for k,v in counts.items():
                    print(f'{indent}{indent}{k} ({v}):',', '.join(sa.printshort() for sa in self.Seqadv if sa.conflict==k))
    def Clone(self,chainmap={},invchainmap={}):
        if len(chainmap)>0:
            newmd=MolData(userMods={},parent_molecule=self.parent_molecule)
            oldmd=self
            for oc in chainmap.keys():
                if oc in oldmd.Chains.keys():
                    pass
                    #print(f'Chain {oc} will be cloned to chain {chainmap[oc]}')
                else:
                    print(f'ERROR: chain {oc} is not in original Chains?')
            # clone Atoms and MissingRes's chain-by-chain
            for oc,nc in chainmap.items():
                newmd.Atoms.extend([a.Clone(chain=nc) for a in oldmd.Atoms if a.chainID==oc])
                newmd.MissingRes.extend([m.Clone(chain=nc) for m in oldmd.MissingRes if m.chainID==oc])
            newmd._MakeResidues()
            # clone SSBonds and XSSBonds
            for oss in oldmd.SSBonds:
                if oss.chainID1 in chainmap and oss.chainID2 in chainmap:
                    nc1=chainmap[oss.chainID1]
                    nc2=chainmap[oss.chainID2]
                    newmd.SSBonds.append(oss.Clone(chainID1=nc1,chainID2=nc2))
            for oxss in oldmd.userXSSBonds:
                if oxss.chainID1 in chainmap and oxss.chainID2 in chainmap:
                    nc1=chainmap[oxss.chainID1]
                    nc2=chainmap[oxss.chainID2]
                    newmd.userXSSBonds.append(oxss.Clone(chainID1=nc1,chainID2=nc2))                    
            # clone links and reprocess
            for olink in oldmd.Links:
                if olink.chainID1 in chainmap and olink.chainID2 in chainmap:
                    nc1=chainmap[olink.chainID1]
                    nc2=chainmap[olink.chainID2]
                    newmd.Links.append(olink.Clone(chainID1=nc1,chainID2=nc2))
            newmd._MakeLinks()
            # copy PDB-explicity sequence discrepancies
            for osa in oldmd.Seqadv:
                if osa.chainID in chainmap:
                    nc=chainmap[osa.chainID]
                    newmd.Seqadv.append(osa.Clone(chain=nc))
            # clone user-specified mutations
            for omut in oldmd.Mutations:
                if omut.chainID in chainmap:
                    nc=chainmap[omut.chainID]
                    newmd.Mutations.append(omut.Clone(chain=nc))
                else:
                    print(f'Error: chainID {omut.chainID} not in chain map?',chainmap)
            #print(f'Clone: {len(newmd.Mutations)} mutations')
            # clone deletions
            for odel in oldmd.Deletions:
                if odel.chainID in chainmap:
                    nc=chainmap[odel.chainID]
                    newmd.Deletions.append(odel.Clone(chain=nc))
            # clone grafts based on target chains
            for ogra in oldmd.Grafts:
                if ogra.target_chain in chainmap:
                    nc=chainmap[ogra.target_chain]
                    newmd.Grafts.append(ogra.Clone(target_chain=nc))
            # clone cleavages based on parent chains
            for oclv in oldmd.Cleavages:
                if oclv.parent_chain in chainmap:
                    nc=chainmap[oclv.parent_chain]
                    newmd.Cleavages.append(oclv.Clone(parent_chain=nc))
            newmd._MakeChains(invchainmap=invchainmap)
            return newmd
        return None

    def GetPatches(self):
        patchstr=''
        SSBonds=self.SSBonds
        XSSBonds=self.userXSSBonds
        nSSBonds=len(SSBonds)-len(XSSBonds)
        patchstr+=f'#### {nSSBonds} SSBONDS:\n'
        for ss in SSBonds:
            if not ss.inlist(XSSBonds):
                patchstr+=ss.psfgen_patchline()
        self._importGrafts()
        Links=self.Links
        patchstr+=f'#### {len(Links)} LINKS:\n'
        for l in Links:
            patchstr+=l.psfgen_patchline()
        return patchstr
    def getGlycanSegnames(self):
        glysegnames=[]
        for s in self.Segments:
            if s.segtype=='GLYGAN':
                glysegnames.append(s.segname)
        return glysegnames
    def MakeSegments(self):
        segments=[]
        for c in self.Chains.values():
            segments.extend(c.MakeSegments(self.Links))
        ''' since this creates segments, we need to update Links so that 
            their patches refer to segments, not chains.  '''
        self._UpdateSegnames()
        self.Segments=segments
        return segments

    def _UpdateSegnames(self):
        for l in self.Links:
            l.updateSegnames(self.Residues)

    def _MakeResidues(self):
        ''' make residues from atoms '''
        self.Residues=[]
        r=0
        for a in self.Atoms:
            if r==0:
                self.Residues.append(Residue(a=a))
                r=self.Residues[-1]
            else:
                if r.resseqnum==a.resseqnum and r.name==a.resname and r.chainID==a.chainID and r.insertion==a.insertion:
                    r.add_atom(a=a)
                else: # begin a new residue
                    self.Residues.append(Residue(a=a))
                    r=self.Residues[-1]
        ''' insert missing residues '''
        for m in self.MissingRes:
            self.Residues.append(Residue(m=m))
    def _MakeChains(self,invchainmap={}):
        self.Chains={}
        for r in self.Residues:
            if r.chainID in self.Chains:
                c=self.Chains[r.chainID]
                c.add_residue(r)
            else:
                if len(invchainmap)>0:
                    source_chainID=invchainmap[r.chainID]
                else:
                    source_chainID=r.chainID
                newChain=Chain(r,parent_molecule=self.parent_molecule,source_chainID=source_chainID)
                self.Chains[newChain.chainID]=newChain
        for c in self.Chains.values():
            c.sort_residues()
#            print(f'Made chain {c.chainID} from source {c.source_chainID}')
        self._ApportionModsToChains()
    def _MakeLinks(self):
        ''' Set all up and down links in residues participating in links '''
        for l in self.Links:
            c1=l.chainID1
            c2=l.chainID2
            #print(l.pdb_line())
            r1=get_residue(self.Residues,c1,l.resseqnum1,l.icode1)
            r2=get_residue(self.Residues,c2,l.resseqnum2,l.icode2)
            l.atom1=get_atom(self.Residues,c1,l.resseqnum1,l.name1,l.icode1)
            l.atom2=get_atom(self.Residues,c2,l.resseqnum2,l.name2,l.icode2)
            l.residue1=r1
            l.residue2=r2
            #print(l.printshort(),str(r1),str(r2))
            # if their chains differ, earlier letters are upstream of later letters
            if _seg_typedict_byresname_[r1.name]=='PROTEIN' and _seg_typedict_byresname_[r2.name]=='GLYCAN':
                r1.linkTo(r2,l)
            elif r1.chainID>r2.chainID:
                r2.linkTo(r1,l)
            elif r1.chainID<r2.chainID:
                r1.linkTo(r2,l)
            else: # they have the same chainID
                if r1.resseqnum>r2.resseqnum:
                    r2.linkTo(r1,l)
                elif r1.resseqnum<r2.resseqnum:
                    r1.linkTo(r2,l)
        for c in self.Chains.values():
            #print('calling group_residues on chain {}'.format(c.chainID))
            c.group_residues()
    def _deleteLink(self,l,deleteDownRes=False):
        #print(f'before deleting link {l.printshort()}: {len(self.Links)} links')
        if l in self.Links:
            r1=l.residue1
            r2=l.residue2
            if r1 in r2.up:
                if deleteDownRes:
         #           print(f'Before recDelDown: {len(self.Residues)} residues, {len(self.Links)} links')
                    _recursiveDelDown(self.Residues,self.Links,r2)
          #          print(f'After recDelDown: {len(self.Residues)} residues, {len(self.Links)} links')
                r1.unlink(r2,l)
            elif r1 in r2.down:
                if deleteDownRes:
           #         print(f'Before recDelDown: {len(self.Residues)} residues, {len(self.Links)} links')
                    _recursiveDelDown(self.Residues,self.Links,r1)
            #        print(f'After recDelDown: {len(self.Residues)} residues, {len(self.Links)} links')
                r2.unlink(r1,l)
            #print(f'Removing link {l.printshort()}')
            self.Links.remove(l)
            #print(f'after deleting link {l.printshort()}: {len(self.Links)} links')
        else:
            print(f'Error: cannot find link {l.printshort()}')
    def _ApportionModsToChains(self):
#        print(f'ApportionMods: {len(self.Mutations)} mutations')
        for c in self.Chains.values():
            cid=c.chainID
            c.Mutations=[m for m in self.Mutations if m.chainID==cid]
            c.Deletions=[d for d in self.Deletions if d.chainID==cid]
            c.Grafts=[g for g in self.Grafts if g.target_chain==cid]
            c.Cleavages=[c for c in self.Cleavages if c.chainID==cid]
    def _checkMutationsAgainstLinks(self):
        mvl=self.userMods['mutationsVsLinks']
        mvs=self.userMods['mutationsVsSSBonds']
        LOC={}
        SSOC={}
        for m in self.Mutations:
            for l in self.Links:
                if m.chainID==l.chainID1 or m.chainID==l.chainID2:
                    if (m.resseqnum==l.resseqnum1 and m.insertion==l.icode1) or (m.resseqnum==l.resseqnum2 and m.insertion==l.icode2):
                        LOC[m]=l
            if m.orig=='CYS':
                for s in self.SSBonds:
                    if m.chainID==s.chainID1 or m.chainID==s.chainID2:
                        if (m.resseqnum==s.resseqnum1 and m.insertion==s.icode1) or (m.resseqnum==s.resseqnum2 and m.insertion==s.icode2):
                            SSOC[m]=s
        reDoChains=False
        if len(LOC)>0:
            #print('#### Mutation-Link conflicts:')
            for k,v in LOC.items():
                #print(k.printshort(),v.printshort(),mvl)
                if mvl=='M':
                    #print('deleting link',v.printshort())
                    self._deleteLink(v,deleteDownRes=True)
                else:
                    self.Mutations.remove(k)
            reDoChains=True
        if len(SSOC)>0:
            #print('#### Mutation-SSBond conflicts:')
            for k,v in SSOC.items():
                #print(k.printshort(),v.printshort(),mvl)
                if mvs=='M':
                    if v in self.SSBonds:  # may have already been removed
                        self.SSBonds.remove(v)
                else:
                    self.Mutations.remove(k)
            reDoChains=True
        if reDoChains:
            self._MakeChains()
#            print('Redid chains: ',', '.join(self.Chains.keys()))
            self.adjustedChainIDs=list(self.Chains.keys())
            self._ApportionModsToChains()
    def _checkDeletionsAgainstLinks(self):
        dvl=self.userMods['deletionsVsLinks']
        dvs=self.userMods['deletionsVsSSBonds']
        LOC={}
        SSOC={}
        for d in self.Deletions:
            for l in self.Links:
                if d.chainID==l.chainID1 or d.chainID==l.chainID2:
                    if (d.resseqnum==l.resseqnum1 and d.insertion==l.icode1) or (d.resseqnum==l.resseqnum2 and d.insertion==l.icode2):
                        LOC[d]=l
            if d.orig=='CYS':
                for s in self.SSBonds:
                    if d.chainID==s.chainID1 or d.chainID==s.chainID2:
                        if (d.resseqnum==s.resseqnum1 and d.insertion==s.icode1) or (d.resseqnum==s.resseqnum2 and d.insertion==s.icode2):
                            SSOC[d]=s
        reDoChains=False
        if len(LOC)>0:
            #print('#### Deletion-Link conflicts:')
            for k,v in LOC.items():
                #print(k.printshort(),v.printshort(),mvl)
                if dvl=='D':
                    #print('deleting link',v.printshort())
                    self._deleteLink(v,deleteDownRes=True)
                else:
                    self.Mutations.remove(k)
            reDoChains=True
        if len(SSOC)>0:
            #print('#### Mutation-SSBond conflicts:')
            for k,v in SSOC.items():
                #print(k.printshort(),v.printshort(),mvl)
                if dvs=='D':
                    if v in self.SSBonds:  # may have already been removed
                        self.SSBonds.remove(v)
                else:
                    self.Mutations.remove(k)
            reDoChains=True
        if reDoChains:
            self._MakeChains()
#            print('Redid chains: ',', '.join(self.Chains.keys()))
            self.adjustedChainIDs=list(self.Chains.keys())
            self._ApportionModsToChains()
    ''' working on below !!! '''
    def _importGrafts(self):
        linksToImport=[]
        linksToEdit=[]
        linksToRemove=[]
        for g in self.Grafts:
            #print('#### importing the following graft into Base as chain {} seg {}'.format(g.ingraft_chainID,g.ingraft_segname,str(g)))
            m=g.molecule
            #for l in m.Links:
            #    l.updateSegnames(m.Residues,m.activeBiologicalAssembly)
            sseg=g.source_segment
            for r in sseg.residues:
                r.segname=g.ingraft_segname
                for a in r.atoms:
                    a.chainID=g.ingraft_chainID
                    a.segname=g.ingraft_segname
                for l in m.Links:
                    if l.isInLink(r.chainID,r.resseqnum):
                        if l.resseqnum1 in g.resid_dict and l.resseqnum2 in g.resid_dict:
                            l.chainID1=g.ingraft_chainID
                            l.chainID2=g.ingraft_chainID
                            l.resseqnum1=g.resid_dict[l.resseqnum1]
                            l.resseqnum2=g.resid_dict[l.resseqnum2]
                            l.segname1=g.ingraft_segname
                            l.segname2=g.ingraft_segname
                            #print('LINK point of import {} {} {} {}'.format(l.segname1,l.resseqnum1,l.segname2,l.resseqnum2))
                            if l not in linksToImport:
                                linksToImport.append(l)
                        #else:
                         #   print('#### this graft link is not imported:')
                          #  print(l.pdb_line())
                         
                r.chainID=g.ingraft_chainID
                r.resseqnum=g.resid_dict[r.resseqnum]
               # print('#### importing graft residue {} into chain {} segment {}'.format(str(r),r.chainID,r.segname))
                self.Residues.append(r)
            # identify Base residues that are grafted over
            graft_overs=[]
            r=get_residue(self.Residues,g.target_chain,g.target_res)
            graft_overs.append(r)
            for rr in r.get_down_group():
                graft_overs.append(rr)
            #print('#### grafting over {:d} base residues rooted at {} in seg {}'.format(len(graft_overs),str(graft_overs[0]),r.segname))
            # find the Base link that joins the target residue to the molecule; target residue is assumed to be the second residue in the link
            for l in self.Links:
                if l.isInLink(graft_overs[0].chainID,graft_overs[0].resseqnum,pos=2):
                    l.chainID2=g.ingraft_chainID
                    l.segname2=g.ingraft_segname
                    l.resseqnum2=g.resid_dict[g.source_res1]
                    linksToEdit.append(l)
            for i,r in enumerate(graft_overs):
               # print('#### graft will remove links involving defunct residue {} in seg {}'.format(str(r),r.segname))
                for l in self.Links:
                    pos=1 if i==0 else ''
                    if l.isInLink(r.chainID,r.resseqnum,pos=pos):
                        if l not in linksToRemove:
                            linksToRemove.append(l)
        if len(linksToImport)>0:
            #print('#### Before import, base has {} links.'.format(len(self.Links)))
            #print('#### The following {} graft links were imported:'.format(len(linksToImport)))
            for l in linksToImport:
                #print(l.pdb_line())
                #print(l.psfgen_patchline())
                self.Links.append(l)
                #print('point of report {} {} {} {}'.format(l.segname1,l.resseqnum1,l.segname2,l.resseqnum2))
            #print('#### Base now has {} links.'.format(len(self.Links)))
            #print('#### The following Base links were edited in-situ')
            for l in linksToEdit:
                #print(l.pdb_line())
                #print(l.psfgen_patchline())
                if l not in self.Links:
                    print('ERROR: In-situ link is not in-situ!')
            #print('#### The following Base links were removed:')
            for l in linksToRemove:
                #print(l.pdb_line())
                #print(l.psfgen_patchline())
                self.Links.remove(l)

