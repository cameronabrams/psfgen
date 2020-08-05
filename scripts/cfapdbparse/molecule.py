from atom import Atom, _PDBAtomNameDict_
from ssbond import SSBond
from missing import Missing
from link import Link
from segment import Segment, _seg_class_
from chain import Chain
from residue import Residue, _PDBResName123_, _pdb_glycans_, _pdb_ions_, _ResNameDict_PDB_to_CHARMM_, _ResNameDict_CHARMM_to_PDB_, get_residue
class Molecule:
    def __init__(self,index,pdb,savepdb=''):
        self.index=index
        self.format='PDB' # default assume this is an actual PDB file from the PDB
        self.RawPDB=[]
        self.keywords=[]
        self.modtyp=[]
        self.title=[]
        self.pdb=pdb
        self.Atoms=[]
        self.Links=[]
        self.SSBonds=[]
        self.MissingRes=[]
        self.MRec={}
        self.Header={}
        self.DBRef={} # outer dictionary referenced by chainID
        self.SeqRes={} # outer: chainID, inner: resnumber, value: resname(PDB)
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
                    #print('i: ',line,end='')
                    #l=self.SSBonds[-1]
                    #print('o: ',l.pdb_line())
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
        self.MakeResidues()
        self.MakeChains()
        self.MakeLinks()
        self.ShowTitle()
        self.ShowHeader()
        self.ShowKeywords()
        self.ShowMasterRecord()
        self.ShowModelType()
        self.ShowDBRef()
        #self.ShowSeqRes()
        if 'CHARMM' in self.keywords:
            self.format=='CHARMM'
            print('### THIS IS A CHARMM-FORMAT PDB FILE')
        print('### Read {:d} pdbrecords from {:s}'.format(len(self.RawPDB),pdb))
    def ParseRemark(self,pdbrecord):
        token=pdbrecord[7:10]
        if token.isdigit():
            code=int(pdbrecord[7:10])
            test_int=pdbrecord[20:26].strip()
            if code==465 and (test_int.isdigit() or (len(test_int)>0 and  test_int[0]=='-')):
                self.MissingRes.append(Missing(pdbrecord))
    def ParseTitle(self,pdbrecord):
        self.title.append(pdbrecord[10:80].strip())
    def ShowTitle(self):
        if len(self.title)>0:
            print('### TITLE records:')
            for i,l in enumerate(self.title):
                print('-> {:1s} {}'.format(' ' if i==0 else str(i),self.title[i]))
        else:
            print('### {} contains no TITLE record'.format(self.pdb))
    def TitleRecord(self):
        retstr=''
        if len(self.title)>0:
           for i,l in enumerate(self.title):
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
        print('### Keywords:')
        for k in self.keywords:
            print('->   {}'.format(k))
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
    def MakeResidues(self):
        self.Residues=[]
        r=0
        for a in self.Atoms:
            if r==0:
                self.Residues.append(Residue(a=a))
                r=self.Residues[-1]
            else:
                if r.resseqnum==a.resseqnum and r.name == a.resname and r.chainID==a.chainID:
                    r.add_atom(a=a)
                else:
                    self.Residues.append(Residue(a=a))
                    r=self.Residues[-1]
        for m in self.MissingRes:
            self.Residues.append(Residue(m=m))
    def MakeChains(self):
        self.Chains=[]
        for r in self.Residues:
            if self.Chains==[]:
                self.Chains.append(Chain(r))
            else:
                found=False
                for c in self.Chains:
                    if c.chainID==r.chainID:
                        found=True
                        c.add_residue(r)
                        break
                if not found:
                    self.Chains.append(Chain(r))
        for c in self.Chains:
            c.sort_residues()
    def MakeLinks(self):
        for l in self.Links:
            r1=get_residue(self.Residues,l.chainID1,l.resseqnum1)
            r2=get_residue(self.Residues,l.chainID2,l.resseqnum2)
            if r1.resseqnum>r2.resseqnum:
                r1.up.append(r2)
                r2.down.append(r1)
            elif r1.resseqnum<r2.resseqnum:
                r1.down.append(r2)
                r2.up.append(r1)
            else:
                if r1.chainID>r2.chainID:
                    r1.up.append(r2)
                    r2.down.append(r1)
                else:
                    r1.down.append(r2)
                    r2.up.append(r1)
        for c in self.Chains:
            c.group_residues() 
    def CleaveChains(self,Cleavages):
        for clv in Cleavages:
            clv_c=-1
            daughter_chainID_ok=True
            for c in self.Chains:
                if c.chainID==clv.parent_chainID:
                    clv_c=c
                if c.chainID==clv.daughter_chainID:
                    daugher_chainID_ok=False
            if clv_c!=-1 and daughter_chainID_ok:
                print('### before cleave:',clv_c)
                daughter=clv_c.Cleave(clv)
                self.Chains.append(daughter)
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
                print('### after cleave:',clv_c,self.Chains[-1])
            else:
                print('### unable to cleave chain {} at position {} to generate {} {}'.format(clv_c.chainID,clv.parent_Cterm_resseqnum,clv.parent_chainID,clv.daughter_chainID))
    def __str__(self):
        return 'Molecule {} {}: {} chains, {} residues, {} atoms, {} links, {} ssbonds'.format(self.index,self.pdb,len(self.Chains),len(self.Residues),len(self.Atoms),len(self.Links),len(self.SSBonds)) 
    def residue_shift(self,chainID,resseqnumshift):
        found=False
        for c in self.Chains:
            if c.chainID==chainID:
                 found=True
                 break;
        if not found:
            print('### Could not apply shift to chain {}: no such chain in Molecule {}'.format(chainID,self.index))
            return -1
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

    def writepsfgeninput(self,fp,Mutations=[],topologies=[]):
        fp.write('if {![info exists PSFGEN_BASEDIR]} {\n'+\
              '    if {[info exists env(PSFGEN_BASEDIR]} {\n'+\
              '        set PSFGEN_BASEDIR $env(PSFGEN_BASEDIR)\n'+\
              '    } else {\n'+\
              '        set PSFGEN_BASEDIR $env(HOME)/research/psfgen\n'+\
              '    }\n'+\
              '}\n')
        fp.write('if {![info exists CHARMM_TOPPARDIR]} {\n'+\
              '    if {[info exists env(CHARMM_TOPPARDIR]} {\n'+\
              '        set TOPPARDIR $env(CHARMM_TOPPARDIR)\n'+\
              '    } else {\n'+\
              '        set TOPPARDIR $env(HOME)/charmm/toppar\n'+\
              '    }\n'+\
              '}\n')
        fp.write('source ${PSFGEN_BASEDIR}/src/loopmc.tcl\n')
        fp.write('source ${PSFGEN_BASEDIR}/scripts/vmdrc.tcl\n')
        fp.write('package require psfgen\n')
        for t in topologies:
            fp.write('topology $TOPPARDIR/{}\n'.format(t))
        fp.write('pdbalias residue HIS HSD\n')
        fp.write('pdbalias atom ILE CD1 CD\n')
        fp.write('pdbalias residue NAG BGNA\n')
        fp.write('pdbalias atom BGNA C7 C\n')
        fp.write('pdbalias atom BGNA O7 O\n')
        fp.write('pdbalias atom BGNA C8 CT\n')
        fp.write('pdbalias atom BGNA N2 N\n')
        fp.write('pdbalias residue SIA ANE5\n')
        fp.write('pdbalias atom ANE5 C10 C\n')
        fp.write('pdbalias atom ANE5 C11 CT\n')
        fp.write('pdbalias atom ANE5 N5 N\n')
        fp.write('pdbalias atom ANE5 O1A O11\n')
        fp.write('pdbalias atom ANE5 O1B O12\n')
        fp.write('pdbalias atom ANE5 O10 O\n')


        for k,v in _ResNameDict_PDB_to_CHARMM_.items():
            fp.write('set RESDICT({}) {}\n'.format(k,v))
        for k,v in _PDBAtomNameDict_.items():
            fp.write('set ANAMEDICT({}) {}\n'.format(k,v))

        fp.write('set logid -1\n')

        fp.write('mol new {}\n'.format(self.pdb))
        Loops=[]
        for c in self.Chains:
            c.MakeSegments(self.Links,Mutations=Mutations)
            for s in c.Segments:
                stan,supp,coor,caco,loops=s.psfgen_segmentstanza()
                fp.write('\n### begin stanza for segment {}\n'.format(s.segname))
                fp.write(supp)
                fp.write(stan+'\n')
                fp.write(coor)
                fp.write(caco)
                if len(loops)>0:
                   Loops.extend(loops)
                fp.write('### end stanza for segment {}\n\n'.format(s.segname))

        fp.write('### SSBONDS:\n')
        for ss in self.SSBonds:
            fp.write(ss.psfgen_patchline())

        fp.write('### LINKS:\n')
        for l in self.Links:
            l.updateSegnames(self.Residues)
            fp.write(l.psfgen_patchline())

        fp.write('guesscoord\n')
        fp.write('regenerate angles dihedrals\n')

        prefix=self.pdb[:]
        prefix=prefix.replace('.pdb','')
        self.psf_outfile='my_{}.psf'.format(prefix)
        self.pdb_outfile='my_{}.pdb'.format(prefix)
        fp.write('writepsf {}\n'.format(self.psf_outfile))
        fp.write('writepdb {}\n'.format(self.pdb_outfile))
        
        # return the list of loops for the post-mod routine to handle
        return Loops

    def Tcl_PrependHeaderToPDB(self,newpdb,psfgen_script):
        hdr='tmp_header.pdb'
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
        # write tcl commands to combine the two files
        psfgen_script.write('exec cat {} >> {}\n'.format(newpdb,hdr))
        psfgen_script.write('exec mv {} {}\n'.format(hdr,newpdb))

