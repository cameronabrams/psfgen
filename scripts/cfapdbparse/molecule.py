from atom import Atom, _PDBAtomNameDict_
from ssbond import SSBond
from missing import Missing
from link import Link
from biomolecule import Biomolecule
from segment import Segment, _seg_class_
from chain import Chain
from seqadv import Seqadv
from mutation import Mutation
from residue import Residue, _PDBResName123_, _pdb_glycans_, _pdb_ions_, _ResNameDict_PDB_to_CHARMM_, _ResNameDict_CHARMM_to_PDB_, get_residue
from revdat import RevDat, FmtDat
from CifFile import ReadCif
_molidcounter_=0
class Molecule:
    def load(self,fp):
        global _molidcounter_
        self.molid=_molidcounter_
        self.molid_varname='m{}'.format(self.molid)
        #print('#### Registered instance of {:s} as molid {:d} with variable name {:s}'.format(self.pdb,self.molid,self.molid_varname))
        self.psfgen_loadstr='mol new {} waitfor all\nset {} [molinfo top get id]\n'.format(self.pdb,self.molid_varname)
        _molidcounter_+=1
        fp.write(self.psfgen_loadstr+'\n')
        
    def __init__(self,pdb='',cif='',isgraft=False,userLinks=[]):
        self.source='RCSB' # default assume this is an actual PDB or CIF file from the RCSB
        self.source_format='CIF' if cif!='' else 'PDB'
        self.molid=-1
        self.molid_varname='UNREGISTERED'
        self.RawPDB=[]
        self.keywords=[]
        self.modtyp=[]
        self.titlelines=[]
        self.Title=''
        self.pdb=pdb
        self.cif=cif
        self.Atoms=[]
        self.Links=userLinks
        self.Chains={} # keyed by chain id 'A', 'B', 'C', etc.
        self.SSBonds=[]
        self.MissingRes=[]
        self.Seqadv=[]
        self.Biomolecules=[]
        self.MRec={}
        self.Header={}
        self.DBRef={} # outer dictionary referenced by chainID
        self.SeqRes={} # outer: chainID, inner: resnumber, value: resname(PDB)
        self.RevDat={}
        if pdb!='':
            with open(pdb) as pdbfile:
                for line in pdbfile:
                    self.RawPDB.append(line)
                    if line[:4] == 'ATOM' or line[:6] == "HETATM":
                        self.Atoms.append(Atom(line))
                    elif line[:4] == 'LINK':
                        self.Links.append(Link(line))
                    elif line[:6] == 'CONECT':
                        #self.Connect.append(Connect(line))
                        pass
                    elif line[:6] == 'SSBOND':
                        self.SSBonds.append(SSBond(line))
                    elif line[:6] == 'SEQADV':
                        self.Seqadv.append(Seqadv(line))
                    elif line[:6] == 'REMARK':
                        self.ParseRemark(line)
                    elif line[:5] == 'TITLE':
                        self.ParseTitle(line)
                    elif line[:6] == 'KEYWDS':
                        self.ParseKeywords(line)
                    elif line[:6] == 'MASTER':
                        self.ParseMasterRecord(line)
                    elif line[:6] == 'HEADER':
                        self.ParseHeader(line)
                    elif line[:6] == 'MDLTYP':
                        self.ParseModelType(line)
                    elif line[:5] == 'DBREF':
                        self.ParseDBRef(line)
                    elif line[:6] == 'SEQRES':
                        self.ParseSeqRes(line)
                    elif line[:6] == 'REVDAT':
                        self.ParseRevisionDate(line)
                    elif line[:6] == 'EXPDTA':
                        self.ParseExpDta(line)
        else:
            print('#### reading {}'.format(cif))
            cf=ReadCif(cif)
            print('#### done')
            db=cf.first_block()
            self.ParseCifDataBlock(db)
        if 'CHARMM' in self.keywords:
            self.source='CHARMM'
        #print('### Read {:d} pdbrecords from {:s}'.format(len(self.RawPDB),pdb))
        self.MakeBiomolecules()
        self.MakeResidues()
        self.MakeChains()
        self.MakeLinks()
        self.MakeSSBonds()
    def summarize(self):
        print('File: {}, Source: {}, Source format: {}'.format(self.pdb if self.source_format=='PDB' else self.cif,self.source,self.source_format))
        print('Title: {}'.format(self.Title))
        self.ShowKeywords()
        if self.source=='RCSB':
            if self.source_format=='PDB':
                print('{}'.format(str(self.FmtDat)))
            else:
                print('CIF Dict. version: {}'.format(self.cif_dict_version))
            print('Last revision: {}'.format(self.ShowRevisions(which='latest',justdates=True)))
            #print('All revisions: {}'.format(self.ShowRevisions(which='all',justdates=False)))
            print('Method: {}; Resolution: {} Ang.'.format(self.ExpDta,self.Resolution))
            print('{} ATOM or HETATOM records.'.format(len(self.Atoms)))
            if len(self.Chains)>0:
               print('Chains: {}'.format(", ".join(c.chainID for c in self.Chains.values())))
            print('Biomolecules:')
            for b in self.Biomolecules:
                b.show()
    def show(self,verbosity):
        print('#'*60)
        print('### MOLID {:s}'.format(self.molid_varname))
        if verbosity>1:
            self.ShowHeader()
            self.ShowKeywords()
            self.ShowMasterRecord()
            self.ShowModelType()
            self.ShowDBRef()
            self.ShowSeqadv(brief=True)
        print('#'*60)
    def ParseRemark(self,pdbrecord):
        token=pdbrecord[7:10].strip()
        if token.isdigit():
            code=int(token)
            if code==465:
                test_int=pdbrecord[20:26].strip()
                if test_int.isdigit() or (len(test_int)>0 and test_int[0]=='-'):
                    self.MissingRes.append(Missing(pdbrecord))
            elif code==350:
                ''' we are reading a BIOMOLECULE stanza '''
                words=pdbrecord.strip().split(' ')
                while words.count(' '):
                    words.remove(' ')
                while words.count(''):
                    words.remove('')
                #print(words)
                if len(words)>2 and words[2].strip()=='BIOMOLECULE:':
                    newBiomolecule=Biomolecule(pdbrecord)
                    self.Biomolecules.append(newBiomolecule)
                elif len(words)>2:
                    ''' do we have a biomolecule instance we are populating? '''
                    if len(self.Biomolecules)>0:
                       bm=self.Biomolecules[-1] # populate newest one
                       bm.parsePDBrecordwords(words)
                #else:
                #    print('ERROR: Non-empty REMARK 350 found outside BIOMOLECULE stanza')
            elif code==2:
                tag=pdbrecord[11:22]
                if tag=='RESOLUTION.':
                    self.Resolution=float(pdbrecord[23:30])
            elif code==4:
                if len(pdbrecord.strip())>11:
                    self.FmtDat=FmtDat(pdbrecord)
                
    def ParseTitle(self,pdbrecord):
        short=pdbrecord[10:80].strip()
        self.titlelines.append(short)
        self.Title=short if len(self.Title)==0 else (self.Title+short)
    def ShowTitleLines(self):
        if len(self.titlelines)>0:
            print('### TITLE records:')
            for i,l in enumerate(self.title):
                print('-> {:1s} {}'.format(' ' if i==0 else str(i),self.titlelines[i]))
        else:
            print('### {} contains no TITLE record'.format(self.pdb))
    def TitleRecord(self):
        retstr=''
        if len(self.titlelines)>0:
           for i,l in enumerate(self.titlelines):
               retstr+='TITLE   {}'.format(' ' if i==0 else str(i))+' {}\n'.format(l)
        return retstr
    def ParseKeywords(self,pdbrecord):
        ctok=pdbrecord[9:10]
        kwdlist=pdbrecord[10:80].split(',')
        for i,k in enumerate(kwdlist):
            if i==0 and ctok.isdigit():
                pk=self.keywords[-1]
                nk=pk+' '+k.strip()
                self.keywords[-1]=nk
            else:
                self.keywords.append(k.strip())
    def ShowKeywords(self):
        print('Keywords: {}'.format(", ".join(self.keywords)))
    def KeywordsRecord(self):
        retstr=''
        if len(self.keywords)>0:
           cpst=', '.join(self.keywords)+' '
           i=0
           sp=[]
           lnlim=69
           while i<len(cpst):
               if cpst[i]==' ':
                   sp.append(i)
               i = i + 1
           nln=len(cpst)//lnlim+1
           beg=0
           for j in range(nln):
               end=-1
               # find rightmost ' ' such that length < lnlim
               for i in range(len(sp)-1,-1,-1):
                   if sp[i]-beg<69:
                      end=sp[i]
                      break
               retstr+='KEYWDS  {:>2s} {}\n'.format(' ' if j==0 else str(j+1),cpst[beg:end])
               beg=end+1
        return retstr
    def ParseModelType(self,pdbrecord):
        ctok=pdbrecord[9:10]
        mdllist=pdbrecord[10:80].split(',')
        for i,m in enumerate(mdllist):
            if i==0 and ctok.isdigit():
                pm=self.modtyp[-1]
                nm=pm+m.strip()
                self.modtyp[-1]=nm
            else:
                self.modtyp.append(m.strip())
    def ShowModelType(self):
        if len(self.modtyp)>0:
            print('### ModelType:')
            for m in self.modtyp:
                print('->   {}'.format(m))
        else:
           print('### {} contains no MODTYP records'.format(self.pdb))
    def ParseMasterRecord(self,pdbrecord):
        self.MRec['numRemark']=int(pdbrecord[10:15])
        self.MRec['numHet']=int(pdbrecord[20:25])
        self.MRec['numHelix']=int(pdbrecord[25:30])
        self.MRec['numSheet']=int(pdbrecord[30:35])
        self.MRec['numSite']=int(pdbrecord[40:45])
        self.MRec['numXform']=int(pdbrecord[45:50])
        self.MRec['numCoord']=int(pdbrecord[50:55])
        self.MRec['numTer']=int(pdbrecord[55:60])
        self.MRec['numConect']=int(pdbrecord[60:65])
        self.MRec['numSeq']=int(pdbrecord[65:70])
    def ShowMasterRecord(self):
        if len(self.MRec)>0:
            print('### MASTER record')
            for k,v in self.MRec.items():
                print('->   {}: {:d}'.format(k,v))
        else:
            print('### {} contains no MASTER record'.format(self.pdb))
    def ParseHeader(self,pdbrecord):
        self.Header['classification']=pdbrecord[10:50]
        self.Header['depDate']=pdbrecord[50:59]
        self.Header['idCode']=pdbrecord[62:66]
    def ShowHeader(self):
        if len(self.Header)>0:
            print('### HEADER record:')
            for k,v in self.Header.items():
                print('->   {}: {}'.format(k,v))
        else:
            print('### {} contains no HEADER record'.format(self.pdb))
    def ParseDBRef(self,pdbrecord):
        if len(pdbrecord)>0:
            chainid=pdbrecord[12:13]
            self.DBRef[chainid]={}
            self.DBRef[chainid]['idCode']=pdbrecord[7:11]
            self.DBRef[chainid]['seqBegin']=int(pdbrecord[14:18])
            self.DBRef[chainid]['insertBegin']=pdbrecord[18:19]
            self.DBRef[chainid]['seqEnd']=int(pdbrecord[20:24])
            self.DBRef[chainid]['insertEnd']=pdbrecord[24:25]
            self.DBRef[chainid]['database']=pdbrecord[26:32]
            self.DBRef[chainid]['dbAccession']=pdbrecord[33:41]
            self.DBRef[chainid]['dbIdCode']=pdbrecord[42:54]
            self.DBRef[chainid]['dbseqBegin']=int(pdbrecord[55:60])
            self.DBRef[chainid]['dbinsBeg']=pdbrecord[60:61]
            self.DBRef[chainid]['dbseqEnd']=int(pdbrecord[62:67])
            self.DBRef[chainid]['dbinsEnd']=pdbrecord[67:69]
    def ShowDBRef(self):
        if len(self.DBRef)>0:
            for ok,ov in self.DBRef.items():
                print('### DBREF {:s}'.format(ok))
                for k,v in ov.items():
                    print('->   {}: {}'.format(k,v))
        else:
            print('### {} contains no DBREF records'.format(self.pdb))
    def ShowBiomolecules(self):
        for b in self.Biomolecules:
            b.show()
    def ParseSeqRes(self,pdbrecord):
        if len(pdbrecord)>0:
            if len(self.DBRef)>0:
                serial=int(pdbrecord[7:10])
                chainid=pdbrecord[11:12]
                if chainid in self.SeqRes:
                    rsn0=serial*13+self.DBRef[chainid]['seqBegin']
                else:
                    self.SeqRes[chainid]={}
                    rsn0=self.DBRef[chainid]['seqBegin']
                rsn=rsn0
                i0=19
                for k in range(13):
                    tok=pdbrecord[i0+k*4:i0+k*4+3].strip()
                    if len(tok)>0:
                        self.SeqRes[chainid][rsn]=tok
                        rsn = rsn + 1
            else:
                print('### Error: {} apparently contains SEQRES records but no DBREF records.'.format(self.pdb))
    def ShowSeqRes(self):
        if len(self.SeqRes)>0:
            for ok,ov in self.SeqRes.items():
                print('### SeqRes for chain {:s}:'.format(ok))
                for k,v in ov.items():
                    print('->   {:d} {:s}'.format(k,v))
        else:
            print('### {} contains no SEQRES records'.format(self.pdb))
    def ParseRevisionDate(self,pdbrecord):
        if len(pdbrecord)>0:
            index=int(pdbrecord[8:10])
            if index in self.RevDat:
                self.RevDat[index].addTokens(pdbrecord)
            else:
                self.RevDat[index]=RevDat(pdbrecord)
    def ShowRevisions(self,which='all',justdates=False):
        maxk=-1
        for k in self.RevDat:
            if k>maxk:
                maxk=k
        if which=='latest' and maxk!=-1:
            if justdates:
                return str(self.RevDat[maxk])
            else:
                return self.RevDat[maxk].show()
        elif which=='all':
            if justdates:
                retstr=''
                for i,r in enumerate(self.RevDat.values()):
                    retstr+='{:s}{}'.format(str(r),', ' if i < len(self.RevDat)-1 else '')
                return retstr 
            else:
                retstr=''
                for r in self.RevDat.values():
                    retstr+=r.show()+'\n'
                return retstr 
           
    def ParseExpDta(self,pdbrecord):
        if len(pdbrecord)>0:
            self.ExpDta=pdbrecord[10:].strip().title()
    def ShowSeqadv(self,brief=False):
        if len(self.Seqadv)>0:
            if not brief:
                for sa in self.Seqadv:
                    print('### SEQADV:',sa)
            else:
                nc=0
                for sa in self.Seqadv:
                    if sa.conflict=='CONFLICT':
                        nc = nc + 1
                print('#### {} contains {:d} SEQADV records including {:d} conflicts'.format(self.pdb,len(self.Seqadv),nc))

    def MakeBiomolecules(self):
        self.chainIDs_allowed=set(['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'])
        chainIDs_detected=set()
        for a in self.Atoms:
            chainIDs_detected.add(a.chainID)
        self.chainIDs_available=sorted(list(self.chainIDs_allowed.difference(chainIDs_detected)))
        chainIDs_detected=sorted(list(chainIDs_detected))
        if len(self.Biomolecules)==0:
            self.Biomolecules.append(Biomolecule('IDENTITY'))
            for c in chainIDs_detected:
                self.Biomolecules[0].chains.append(c)
            #print('#### added identity biomt to pdb without any',self.Biomolecules[0].chains)
        for b in self.Biomolecules:
            for c in b.chains:
                for t in b.biomt:
                    new_chainID=c
                    if not t.isidentity():
                        new_chainID=self.chainIDs_available.pop(0)
                    t.mapchains(c,new_chainID)
    def GetBiomoleculeByChain(self,c):
        for b in self.Biomolecules:
            if c in b.chains:
                return b
        return None
    def MakeResidues(self):
        ''' make residues from atoms '''
        self.Residues=[]
        r=0
        for a in self.Atoms:
            if r==0:
                self.Residues.append(Residue(a=a))
                r=self.Residues[-1]
            else:
                if r.resseqnum==a.resseqnum and r.name == a.resname and r.chainID==a.chainID:
                    r.add_atom(a=a)
                else: # begin a new residue
                    self.Residues.append(Residue(a=a))
                    r=self.Residues[-1]
        ''' insert missing residues '''
        for m in self.MissingRes:
            self.Residues.append(Residue(m=m))
    def MakeChains(self):
        for r in self.Residues:
            if r.chainID in self.Chains:
                c=self.Chains[r.chainID]
                c.add_residue(r)
            else:
                newChain=Chain(r,parent_molecule=self)
                self.Chains[newChain.chainID]=newChain
        for c in self.Chains.values():
            c.sort_residues()
            #print('#### chain {} parent_molecule {}'.format(c.chainID,c.parent_molecule.molid_varname))
        #print('### Created {} chains.'.format(len(self.Chains)))
    def MakeSegments(self,userMutations=[],userGrafts=[],userAttachments=[]):
        nseg=0
        for c in self.Chains.values():
             tmut=[_ for _ in userMutations if _.chainID==c.chainID]
             tgra=[_ for _ in userGrafts if _.target_chain==c.chainID]
             tatt=[_ for _ in userAttachments if _.target_chain==c.chainID]
             c.MakeSegments(self.Links,tmut,tgra,tatt)
             nseg+=len(c.Segments)
        print('### Created {} segments.'.format(nseg))
    def MakeSSBonds(self):
        #print('#### MakeSSBonds is working with {} bonds'.format(len(self.SSBonds)))
        ''' we need to replicate all SSBOND across new chains in the list of Biomolecules '''
        if len(self.Biomolecules)>0:
            newssbonds=[]
            for b in self.Biomolecules:
                for t in b.biomt:
                    if not t.isidentity():
                        #print('#### we are replicating SSBONDS!')
                        for ss in self.SSBonds:
                            newc1=t.get_replica_chainID(ss.chainID1)
                            newc2=t.get_replica_chainID(ss.chainID2)
                            if newc1!=ss.chainID1 or newc2!=ss.chainID2:
                                bb=SSBond(ss.pdbrecord)
                                bb.chainID1=newc1
                                bb.chainID2=newc2
                                newssbonds.append(bb)
            self.SSBonds.extend(newssbonds)
       # print('#### MakeSSBonds finishes with {} bonds.'.format(len(self.SSBonds)))
 
    def MakeLinks(self):
        ''' Set all up and down links in residues participating in links '''
        for l in self.Links:
            r1=get_residue(self.Residues,l.chainID1,l.resseqnum1)
            r2=get_residue(self.Residues,l.chainID2,l.resseqnum2)
            # if their chains differ, earlier letters are upstream of later letters
            if _seg_class_[r1.name]=='PROTEIN' and _seg_class_[r2.name]=='GLYCAN':
                r1.down.append(r2)
                r2.up.append(r1)
            elif r1.chainID>r2.chainID:
                r1.up.append(r2)
                r2.down.append(r1)
                #print(r2,'->',r1)
            elif r1.chainID<r2.chainID:
                r1.down.append(r2)
                r2.up.append(r1)
                #print(r1,'->',r2)
            else: # they have the same chainID
               if r1.resseqnum>r2.resseqnum:
                   r1.up.append(r2)
                   r2.down.append(r1)
                   #print(r2,'->',r1)
               elif r1.resseqnum<r2.resseqnum:
                   r1.down.append(r2)
                   r2.up.append(r1)
                   #print(r1,'->',r2)
        for c in self.Chains.values():
            #print('calling group_residues on chain {}'.format(c.chainID))
            c.group_residues()
        ''' make additional link copies if there are biomt transformations '''
        if len(self.Biomolecules)>0:
            newlinks=[]
            for b in self.Biomolecules:
                for t in b.biomt:
                    if not t.isidentity():
                        for l in self.Links:
                            newc1=t.get_replica_chainID(l.chainID1)
                            newc2=t.get_replica_chainID(l.chainID2)
                            if newc1!=l.chainID1 or newc2!=l.chainID2:
                                ll=Link(l.pdbrecord)
                                ll.chainID1=newc1
                                ll.chainID2=newc2
                                newlinks.append(ll)
            self.Links.extend(newlinks)
    def CleaveChains(self,Cleavages):
        if len(self.chainIDs_available)<len(Cleavages):
            print("### WARNING: insufficient chainID's are available for {} cleavages".format(len(Cleavages)))
            print("### {} chain ID's available:".format(len(self.chainIDs_available)),self.chainIDs_available)
            return -1
        for clv in Cleavages:
            daughter_chain_ok=False
            if clv.parent_chainID in self.Chains:
                clv_c=self.Chains[clv.parent_chainID]
                clv.daughter_chainID=self.chainIDs_available.pop(0)
                daughter=clv_c.Cleave(clv)
                self.Chains[daughter.chainID]=daughter
                b=self.GetBiomoleculeByChain(clv.parent_chainID)
                b.chains.append(daughter.chainID)
                for s in self.SSBonds:
                    if s.chainID1==clv_c.chainID and s.resseqnum1>clv.parent_Cterm_resseqnum:
                        s.chainID1=daughter.chainID
                    if s.chainID2==clv_c.chainID and s.resseqnum2>clv.parent_Cterm_resseqnum:
                        s.chainID2=daughter.chainID
                for l in self.Links:
                    if l.chainID1==clv_c.chainID and daughter.has_resseqnum(l.resseqnum1):
                        l.chainID1=daughter.chainID
                    if l.chainID2==clv_c.chainID and daughter.has_resseqnum(l.resseqnum2):
                        l.chainID2=daughter.chainID
            else:
                print('### unable to cleave chain {} at position {}'.format(clv_c.chainID,clv.parent_Cterm_resseqnum))
    def __str__(self):
        return 'Molecule {}: {} chains, {} residues, {} atoms, {} links, {} ssbonds'.format(self.pdb,len(self.Chains),len(self.Residues),len(self.Atoms),len(self.Links),len(self.SSBonds)) 
    def residue_shift(self,chainID,resseqnumshift):
        if chainID not in self.Chains:
            print('#### Warning: cannot shift in chain {}: no such chain'.format(chainID))
            return -1
        c=self.Chains[chainID]
        for r in c.residues:
            r.residue_shift(resseqnumshift)
        for l in self.Links.L:
            if l.chainID1==chainID:
                l.resseqnum1+=resseqnumshift
            if l.chainID2==chainID:
                l.resseqnum2+=resseqnumshift
        for ss in self.SSBonds:
            if ss.chainID1==chainID:
                ss.resseqnum1+=resseqnumshift
            if ss.chainID2==chainID:
                ss.resseqnum2+=resseqnumshift
        return 0

    def WritePsfgenInput(self,fp,userMutations=[],prefix='my_',fixConflicts=False,fixEngineeredMutations=False,userGrafts=[],userAttach=[],userSSBonds=[],userIgnoreChains=[],includeTerminalLoops=False,removePDBs=True):

        self.load(fp)
        for g in userGrafts:
            g.load(fp)
        for a in userAttach:
            pass # a.load(fp)

        fp.write('mol top ${}\n'.format(self.molid_varname))

        ''' update the user-specified list of mutations:
            1. include any conflict-fixes specified in the RCSB-format PDB file
            2. include any possible replicas created when new chains are eventually 
               created via BIOMT operations
        '''
        newmutations=[]
        for sa in self.Seqadv:
            if fixConflicts==True and sa.conflict=='CONFLICT':
                newmutations.append(Mutation(seqadv=sa))
            elif fixEngineeredMutations==True and sa.conflict=='ENGINEERED MUTATION':
                newmutations.append(Mutation(seqadv=sa))
        userMutations.extend(newmutations)
        replica_mutations=[]
        for b in self.Biomolecules:
            for t in b.biomt:
                for m in userMutations:
                    newc=t.get_replica_chainID(m.chainID)
                    if newc!=m.chainID:
                        mm=m.replicate(newchainID=t.get_replica_chainID(m.chainID))
                        replica_mutations.append(mm)
        userMutations.extend(replica_mutations)

        Loops=[]
        for c in self.Chains.values():
            if c.chainID not in userIgnoreChains:
                #print('#### segmentifying chain {}'.format(c.chainID))
                b=self.GetBiomoleculeByChain(c.chainID)
                #print('#### Chain {} is claimed by Biomolecule {}'.format(c.chainID,b.index))
                c.MakeSegments(self.Links,Mutations=userMutations,Grafts=userGrafts,Attachments=userAttach)
                for s in c.Segments:
                    for t in b.biomt:
                        #print('#### Chain {} replica in tmat {:d}: {}'.format(c.chainID,t.index,t.get_replica_chainID(c.chainID)))
                        stanza,loops=s.write_psfgen_stanza(includeTerminalLoops=includeTerminalLoops,tmat=t)
                        Loops.extend(loops) # psfgen postprocessing needs loop info
                        fp.write('\n#### Begin stanza for segment {} biomolecule {} tmat {}\n'.format(s.segname,b.index,t.index))
                        fp.write(stanza)
                        fp.write('#### End stanza for segment {}\n\n'.format(s.segname))
                    for p in s.pdbfiles:
                        if removePDBs:
                            fp.write('file delete {}\n'.format(p))
        fp.write('#### SSBONDS:\n')
        for ss in self.SSBonds:
            if ss.chainID1 not in userIgnoreChains and ss.chainID2 not in userIgnoreChains:
                fp.write(ss.psfgen_patchline())
        for ss in userSSBonds:
            fp.write(ss.psfgen_patchline())

        if len(userGrafts)>0:
            print('#### Importing {} grafts'.format(len(userGrafts)))
            self.importGrafts(userGrafts)
        fp.write('#### {} LINKS:\n'.format(len(self.Links)))
        for l in self.Links:
            ''' The generation of new segments for glycans and/or cleavage necessitates
                updating segment names in the LINK records read from the original PDB.
                updateSegnames() reassigns residue1/2, atom1/2, and segname1/2 
                attributes of each link instance. '''
            l.updateSegnames(self.Residues,self.Biomolecules)
            fp.write(l.psfgen_patchline())

        for b in self.Biomolecules:
            pass

        fp.write('#### END OF SEGMENTS\n')

        fp.write('guesscoord\n')
        fp.write('regenerate angles dihedrals\n')

        code=self.pdb[:]
        code=code.replace('.pdb','')
        self.psf_outfile='{}{}.psf'.format(prefix,code)
        self.pdb_outfile='{}{}.pdb'.format(prefix,code)
        fp.write('writepsf {}\n'.format(self.psf_outfile))
        fp.write('writepdb {}\n'.format(self.pdb_outfile))
        
        # return the list of loops for the post-mod routine to handle
        return Loops

    def importGrafts(self,userGrafts):
        linksToImport=[]
        linksToEdit=[]
        linksToRemove=[]
        for g in userGrafts:
            #print('#### importing the following graft into Base as chain {} seg {}'.format(g.ingraft_chainID,g.ingraft_segname,str(g)))
            m=g.molecule
            for l in m.Links:
                l.updateSegnames(m.Residues,m.Biomolecules)
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
    def Tcl_PrependHeaderToPDB(self,newpdb,psfgen_script,hdr='charmm_header.pdb'):
        _tmpfile_='_tmpfile_'
        fp=open(hdr,'w')
        fp.write(self.TitleRecord())
        self.keywords.append('CHARMM')
        fp.write(self.KeywordsRecord())
        #write ssbonds and links
        for ss in self.SSBonds:
            fp.write(ss.pdb_line()+'\n')
        for l in self.Links:
            fp.write(l.pdb_line()+'\n')
        fp.close()
        # write tcl commands to combine the two files, preserving the header for other uses
        psfgen_script.write('exec cat {} {} > {}\n'.format(hdr,newpdb,_tmpfile_))
        psfgen_script.write('exec mv {} {}\n'.format(_tmpfile_,newpdb))
        psfgen_script.write('exec rm -f {}\n'.format(_tmpfile_))

    def ParseCifDataBlock(self,db):
        structs={}
        for k in db.keys():
            kk=k.split('.')
            s=kk[0]
            if s in structs:
                i+=1
                structs[s][kk[1]]=i
            else:
                i=0
                structs[s]={}
                structs[s][kk[1]]=i
        self.Title=db['_struct.title']
        self.keywords=[_.strip() for _ in db['_struct_keywords.text'].split(',')]
        self.cif_dict_version=db['_audit_conform.dict_version']
        self.ExpDta=db['_exptl.method'].title()
        self.Resolution=db['_refine.ls_d_res_high']
        x=db.GetLoop('_pdbx_audit_revision_history.ordinal')
        keys=[_ for _ in structs['_pdbx_audit_revision_history'].keys()]
        for y in x:
            d={}
            for v in keys:
                d[v]=y[structs['_pdbx_audit_revision_history'][v]]
            self.RevDat[int(d['ordinal'])]=RevDat(d,fmt='CIF')
        isloop=True if db.FindLoop('_pdbx_struct_assembly.id')!=-1 else False
        if isloop:
            pass
        else:
            self.Biomolecules.append(Biomolecule(cifdb=db))


