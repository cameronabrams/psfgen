from atom import Atom, _PDBAtomNameDict_
from ssbond import SSBond
from missing import Missing
from link import Link
from biomolecule import Biomolecule
from segment import Segment, _seg_typedict_byresname_
from chain import Chain
from seqadv import Seqadv
from mutation import Mutation
from residue import Residue, _PDBResName123_, _pdb_glycans_, _pdb_ions_, _ResNameDict_PDB_to_CHARMM_, _ResNameDict_CHARMM_to_PDB_, get_residue
from revdat import RevDat, FmtDat
from CifFile import ReadCif
from cifutil import *
from moldata import MolData
_molidcounter_=0
class Molecule:
    def load(self,fp):
        global _molidcounter_
        self.molid=_molidcounter_
        self.molid_varname='m{}'.format(self.molid)
        #print('#### Registered instance of {:s} as molid {:d} with variable name {:s}'.format(self.pdb,self.molid,self.molid_varname))
        self.psfgen_loadstr='mol new {} waitfor all\nset {} [molinfo top get id]\n'.format(self.pdb,self.molid_varname)
        if self.cif!=None:
            # need to renumber and rechain to user specs
            self.psfgen_loadstr+='set ciftmp [atomselect ${} all]\n'.format(self.molid_varname)
            self.psfgen_loadstr+='$ciftmp set chain [list {}]\n'.format(" ".join([_.chainID for _ in self.Atoms]))
            self.psfgen_loadstr+='$ciftmp set resid [list {}]\n'.format(" ".join([str(_.resseqnum) for _ in self.Atoms]))
        _molidcounter_+=1
        fp.write(self.psfgen_loadstr+'\n')
        return self.molid
    def ReadPDB(self,pdb):
        with open(pdb) as pdbfile:
            for line in pdbfile:
                self.RawPDB.append(line)
                key=line[:6].strip()
                if key in self.md.PDBClonables:
                    self.md.Parse(line)
                elif key == 'REMARK':
                    if line[7:10]=='465':
                        self.md.Parse(line)
                    else:
                        self.ParseRemark(line)
                elif key == 'TITLE':
                    self.ParseTitle(line)
                elif key == 'KEYWDS':
                    self.ParseKeywords(line)
                elif key == 'MASTER':
                    self.ParseMasterRecord(line)
                elif key == 'HEADER':
                    self.ParseHeader(line)
                elif key == 'MDLTYP':
                    self.ParseModelType(line)
                elif key == 'DBREF':
                    self.ParseDBRef(line)
                elif key == 'SEQRES':
                    self.ParseSeqRes(line)
                elif key == 'REVDAT':
                    self.ParseRevisionDate(line)
                elif key == 'EXPDTA':
                    self.ParseExpDta(line)
                else:
                    self.unusedKeys.add(key)

    def __init__(self,pdb=None,cif=None,userMods={},requestedBiologicalAssembly=None):
        self.source='RCSB' # default assume this is an actual PDB or CIF file from the RCSB
        self.source_format='CIF' if cif!=None else 'PDB'
        self.molid=-1
        self.molid_varname='UNREGISTERED'
        self.md=MolData(userMods=userMods,parent_molecule=self)
        self.RawPDB=[]
        self.keywords=[]
        self.modtyp=[]
        self.titlelines=[]
        self.Title=''
        self.pdb=pdb if pdb!=None else cif
        self.unusedKeys=set()
        self.cif=cif
        self.Biomolecules=[]
        self.activeBiologicalAssembly=None
        self.MRec={}
        self.Header={}
        self.DBRef={} # outer dictionary referenced by chainID
        self.SeqRes={} # outer: chainID, inner: resnumber, value: resname(PDB)
        self.RevDat={}
        ltrs=set()
        for x in range(26):
            ltrs.add(chr(ord('A')+x))
            ltrs.add(chr(ord('a')+x))
        self.chainIDs_allowed=ltrs
        if pdb!=None:
            self.ReadPDB(pdb)
        elif cif!=None:
            cf=ReadCif(cif)
            db=cf.first_block()
            self.ParseCifDataBlock(db)
        else:
            print('Error: Molecule __init__ called without a record.')
            return
        if 'CHARMM' in self.keywords:
            self.source='CHARMM'
        self.md.MakeConstructs()
        self.MakeBiomolecules(requestedBiologicalAssembly=requestedBiologicalAssembly)
    def summarize(self):
        print('File: {}, Source: {}, Source format: {}'.format(self.pdb if self.source_format=='PDB' else self.cif,self.source,self.source_format))
        print('   Title: {}'.format(self.Title))
        self.ShowKeywords(indent=' '*3)
        if self.source=='RCSB':
            if self.source_format=='PDB':
                print('   {}'.format(str(self.FmtDat)))
            else:
                print('   CIF Dict. version: {}'.format(self.cif_dict_version))
            print('   Last revision: {}'.format(self.ShowRevisions(which='latest',justdates=True)))
            #print('All revisions: {}'.format(self.ShowRevisions(which='all',justdates=False)))
            print('   Method: {}; Resolution: {} Ang.'.format(self.ExpDta,self.Resolution))
            print('   {} ATOM or HETATOM records.'.format(len(self.md.Atoms)))
            print('   {} unique residues, {} missing.'.format(len(self.md.Residues),len(self.md.MissingRes)))
            print('   {} disulfides; {} covalent links.'.format(len(self.md.SSBonds),len(self.md.Links)))
            self.md.ShowSeqadv(brief=True)
            if len(self.md.Chains)>0:
               print('   {} chains: {}'.format(len(self.md.Chains),", ".join(c.chainID for c in self.md.Chains.values())))
            print('   {} Biological assemblies:'.format(len(self.Biomolecules)))
            for b in self.Biomolecules:
                b.show(indent=' '*6,isActive=(b is self.activeBiologicalAssembly))
            print('')
    def show(self,verbosity):
        print('#'*60)
        print('### MOLID {:s}'.format(self.molid_varname))
        if verbosity>1:
            self.ShowHeader()
            self.ShowKeywords()
            self.ShowMasterRecord()
            self.ShowModelType()
            self.ShowDBRef()
            self.md.ShowSeqadv(brief=True)
        print('#'*60)
    def ParseRemark(self,pdbrecord):
        token=pdbrecord[7:10].strip()
        if token.isdigit():
            code=int(token)
            if code==465:
                ''' we are reading a missing residue stanza '''
                test_int=pdbrecord[20:26].strip()
                if test_int.isdigit() or (len(test_int)>0 and test_int[0]=='-'):
                    self.md.MissingRes.append(Missing(pdbrecord))
            elif code==350:
                ''' we are reading a BIOMOLECULE stanza '''
                words=pdbrecord.strip().split(' ')
                while words.count(' '):
                    words.remove(' ')
                while words.count(''):
                    words.remove('')
                #print(words)
                if len(words)>2 and words[2].strip()=='BIOMOLECULE:':
                    newBiomolecule=Biomolecule(pdbrecord,parent_molecule=self)
                    self.Biomolecules.append(newBiomolecule)
                elif len(words)>2:
                    ''' do we have a biomolecule instance we are populating? '''
                    if len(self.Biomolecules)>0:
                       bm=self.Biomolecules[-1] # populate newest one
                       #print('parsing biomolecule using ',words)
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
    def ShowKeywords(self,indent=''):
        print('{}Keywords: {}'.format(indent,', '.join(self.keywords)))
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

    def MakeBiomolecules(self,requestedBiologicalAssembly=None):
        self.Biomolecules.insert(0,Biomolecule(parent_molecule=self)) # make the 0th biomolecule the asymmetric unit
        asymmetricUnit=self.Biomolecules[0]
        asymmetricUnit.initializeAsymmetricUnit(self.md)
        self.ChainIDDepot=sorted(list(self.chainIDs_allowed.difference(set(asymmetricUnit.chainIDs))))
        if requestedBiologicalAssembly==None:
            self.activeBiologicalAssembly=self.Biomolecules[0]
        else:
            nBA=len(self.Biomolecules)
            if requestedBiologicalAssembly>nBA:
                print('Error: requested biological assembly index {} exceeds number of available assemblies {}.'.format(requestedBiologicalAssembly,nBA))
                print('Using the asymmetric unit.')
                self.activeBiologicalAssembly=self.Biomolecules[0]
            else:
                self.activeBiologicalAssembly=self.Biomolecules[requestedBiologicalAssembly]
                b=self.activeBiologicalAssembly
                if not b is asymmetricUnit:
                    self.ChainIDDepot=b.inheritConstructs(asymmetricUnit,self.ChainIDDepot)
    def getGlycanSegnames(self):
        glysegnames=[]
        for t in self.activeBiologicalAssembly.biomt:
            glysegnames.extend(t.md.getGlycanSegnames())
        return glysegnames

    def CleaveChains(self,Cleavages):
        if len(self.chainIDs_available)<len(Cleavages):
            print("### WARNING: insufficient chainID's are available for {} cleavages".format(len(Cleavages)))
            print("### {} chain ID's available:".format(len(self.chainIDs_available)),self.chainIDs_available)
            return -1
        for clv in Cleavages:
            # daughter_chain_ok=False
            if clv.parent_chainID in self.activeBiologicalAssembly.activeChainIDs:
                clv_c=self.Chains[clv.parent_chainID]
                clv.daughter_chainID=self.chainIDs_available.pop(0)
                daughter=clv_c.Cleave(clv)
                self.Chains[daughter.chainID]=daughter
                self.activeBiologicalAssembly.activeChainIDs.append(daughter.chainID)
                b=self.activeBiologicalAssembly.getBiomT(chainID=clv.parent_chainID)
                b.chainIDs.append(daughter.chainID)
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
        return 'Molecule {}: {} chains, {} residues, {} atoms, {} links, {} ssbonds, {} biological assemblies'.format(self.pdb,len(self.md.Chains),len(self.md.Residues),len(self.md.Atoms),len(self.md.Links),len(self.md.SSBonds),len(self.Biomolecules)) 
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

    def WritePsfgenInput(self,fp,prefix='DEFAULTPREFIX'):
        ''' issue the mol load commands '''
        molids=[self.load(fp)]
        for t in self.activeBiologicalAssembly.biomt:
            for g in t.md.Grafts:
                molids.append[g.load(fp)]
        print('Molids:',molids)
        fp.write(f'mol top ${self.molid_varname}\n')
        fp.write('#### BEGIN SEGMENTS\n')
        Loops=[]
        for t in self.activeBiologicalAssembly.biomt:
            Segments=t.md.MakeSegments()
            for s in Segments:
                stanza,loops=s.psfgen_stanza(includeTerminalLoops=self.md.userMods['includeTerminalLoops'],tmat=t)
                Loops.extend(loops) # psfgen postprocessing needs loop info
                fp.write(stanza)
        fp.write('#### END SEGMENTS\n')
        fp.write('#### BEGIN PATCHES\n')
        for t in self.activeBiologicalAssembly.biomt:
            fp.write(t.md.GetPatches())
        fp.write('#### END PATCHES\n')

        fp.write('guesscoord\n')
        fp.write('regenerate angles dihedrals\n')

        code='.'.join(self.pdb.split('.')[:-1])
#        code=self.pdb[:]
#        code=code.replace('.pdb','')
        self.psf_outfile='{}{}.psf'.format(prefix,code)
        self.pdb_outfile='{}{}.pdb'.format(prefix,code)
        fp.write('writepsf cmap {}\n'.format(self.psf_outfile))
        fp.write('writepdb {}\n'.format(self.pdb_outfile))
        
        # return the list of loops for the post-mod routine to handle
        return Loops

    def Tcl_PrependHeaderToPDB(self,newpdb,psfgen_script,hdr='charmm_header.pdb'):
        _tmpfile_='_tmpfile_'
        fp=open(hdr,'w')
        fp.write(self.TitleRecord())
        self.keywords.append('CHARMM')
        fp.write(self.KeywordsRecord())
        #write ssbonds and links
        for t in self.activeBiologicalAssembly.biomt:
            for ss in t.md.SSBonds:
                fp.write(ss.pdb_line()+'\n')
            for l in t.md.Links:
                fp.write(l.pdb_line()+'\n')
        fp.close()
        # write tcl commands to combine the two files, preserving the header for other uses
        psfgen_script.write('exec cat {} {} > {}\n'.format(hdr,newpdb,_tmpfile_))
        psfgen_script.write('exec mv {} {}\n'.format(_tmpfile_,newpdb))
        psfgen_script.write('exec rm -f {}\n'.format(_tmpfile_))

    def ParseCifDataBlock(self,db):
        structs=CIFMakeStructs(db)
        self.Title=db['_struct.title']
        self.keywords=[_.strip() for _ in db['_struct_keywords.text'].split(',')]
        self.cif_dict_version=db['_audit_conform.dict_version']
        self.ExpDta=db['_exptl.method'].title()
        self.Resolution=db['_refine.ls_d_res_high']
        
        struk='_pdbx_audit_revision_history'
        dlist=GetCIFStructDictList(db,structs,struk)
        for d in dlist:
            self.RevDat[int(d['ordinal'])]=RevDat(d,fmt='CIF')
        
        # "assemblies" in the CIF are "biomolecules" in the PDB
        struk='_pdbx_struct_assembly'
        dlist=GetCIFStructDictList(db,structs,struk)
        for d in dlist:
            self.Biomolecules.append(Biomolecule(cifdict=d))
        self.CIFParseBiomolecules(GetCIFStructDictList(db,structs,'_pdbx_struct_assembly_gen'),GetCIFStructDictList(db,structs,'_pdbx_struct_oper_list'))
        dlist=GetCIFStructDictList(db,structs,'_atom_site')
        self.CIFParseAtoms(dlist)
        dlist=GetCIFStructDictList(db,structs,'_pdbx_unobs_or_zero_occ_residues')
        self.CIFParseMissing(dlist)
        dlist=GetCIFStructDictList(db,structs,'_struct_conn')
        self.CIFParseConnections(dlist)

    def CIFParseBiomolecules(self,genl,operl):
        # for each gen, associate with a biomoleculr  
        #print(genl)
        #print(operl)
        for g in genl:
            index=int(g['assembly_id'])-1
            operids=g['oper_expression'].split(',')
            chains=g['asym_id_list'].split(',')
            for i in operids:
                useme={}
                for od in operl:
                    if od['id']==i:
                        useme=od
                if len(useme)==0:
                    print('Error: oper index {} not found in oper_list'.format(i))
                else:
                    self.Biomolecules[index].CIFBiomT(useme)
                    self.Biomolecules[index].chains=chains[:]

    def CIFParseAtoms(self,alist):
        for a in alist:
            self.Atoms.append(Atom(cifdict=a))

    def CIFParseMissing(self,mlist):
        for m in mlist:
            self.MissingRes.append(Missing(cifdict=m))

    def CIFParseConnections(self,clist):
        for c in clist:
            typ=c['conn_type_id']
            if typ=='disulf':
                self.SSBonds.append(SSBond(cifdict=c))
            elif typ=='covale':
                self.Links.append(Link(cifdict=c))


