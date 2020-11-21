class Biomolecule:
    ''' Container for handling info for "REMARK 350 BIOMOLECULE: #" stanzas in RCSB PDB files
        or _pdbx_struct blocks in RCSB CIF files '''
    def __init__(self,pdbrecord=None,cifdict=None):
        #print('__init__ with {}'.format(pdbrecord))
        self.chains=[]
        self.biomt=[]
        self.pdbx_struct={}
        if pdbrecord!=None and cifdict!=None:
           print('Warning: Biomolecule.__init__ called with both a pdbrecord and cifdb; using pdbrecord')
        if pdbrecord!=None:
            if 'BIOMOLECULE:' in pdbrecord:
               self.index=int(pdbrecord[23:25].strip())
            if 'IDENTITY' in pdbrecord:  
                # caller would like an identity created, most likely because no 
                # BIOMOLECULE or PDBX_STRUCT was in the input file
                self.index=1
                self.biomt.append(BiomT())
        elif cifdict!=None:
            self.index=int(cifdict['id'])
            for k in ['method_details','oligomeric_details','oligomeric_count']:
                self.pdbx_struct[k]=cifdict[k]
            
    def parsePDBrecordwords(self,words):
        if len(words)>2:
            phrase=' '.join(words[2:])
          #  print(phrase)
            if 'AUTHOR DETERMINED BIOLOGICAL UNIT:' in phrase:
           #     print(words[-1])
                self.pdbx_struct['author_biol_unit']=words[-1]
            elif 'SOFTWARE DETERMINED QUATERNARY STRUCTURE:' in phrase:
                self.pdbx_struct['software_quat_struct']=words[-1]
            elif 'SOFTWARE USED:' in phrase:
                self.pdbx_struct['software_used']=[_.replace(',','') for _ in words[4:]]
            elif 'TOTAL BURIED SURFACE AREA:' in phrase:
                self.pdbx_struct['total_buried_surface_area']=int(words[-2])
                self.pdbx_struct['surface_area_units']=words[-1]
            elif 'SURFACE AREA OF THE COMPLEX' in phrase:
                self.pdbx_struct['surface_area_complex']=int(words[-2])
            elif 'CHANGE IN SOLVENT FREE ENERGY:' in phrase:
                self.pdbx_struct['change_solv_fe']=float(words[-2])
                self.pdbx_struct['fe_units']=words[-1]
            elif 'APPLY THE FOLLOWING TO CHAINS:' in phrase:
                self.chains.extend([_.replace(',','') for _ in words[7:]])
            elif 'AND CHAINS:' in phrase:
                self.chains.extend([_.replace(',','') for _ in words[4:]])
            elif 'BIOMT' in words[2]:
                self.parseBIOMT(words)
            else:
                print('#### ERROR: Unrecognized PDB-format REMARK 350 line:')
                print(' '.join(words))
    def show(self,indent='    '):
        print('{}Biomolecule {:d}'.format(indent,self.index))
        print(self.pdbx_struct,self.chains)
        for b in self.biomt:
            b.show(indent)
    def parseBIOMT(self,words):
        ax=int(words[2][-1])
        if ax==1:
            self.biomt.append(BiomT())
        self.biomt[-1].parseBIOMT(ax,words)
    def CIFBiomT(self,cifdict):
        self.biomt.append(BiomT())
        self.biomt[-1].CIFBiomT(cifdict)

class BiomT:
    def __init__(self):
        self.index=-1
        self.tmat=[[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 1, 0]]
        self.replicachainID_from_sourcechainID={}
        self.sourcechainID_from_replicachainID={}
    def parseBIOMT(self,ax,words):
        if self.index==-1:
            self.index=int(words[3])
        vals=[]
        for w in words[4:]:
           vals.append(float(w))
        self.tmat[ax-1]=vals
    def CIFBiomT(self,cifdict):
        self.index=int(cifdict['id'])
        for i in range(3):
            for j in range(3):
                self.tmat[i][j]=float(cifdict['matrix[{}][{}]'.format(i+1,j+1)])
            self.tmat[i][3]=float(cifdict['vector[{}]'.format(i+1)])

    def show(self,indent='    '):
        print('{}BIOMT {:d}'.format(indent,self.index))
        if not self.isidentity()
            print('{}    TMAT'.format(indent),self.tmat)
            print('{}    REPC'.format(indent),self.replicachainID_from_sourcechainID)
        else:
            print('{}    IDENTITY'.format(indent))
    def isidentity(self):
        t=self.tmat
        if t[0][0]==1.0 and t[1][1]==1.0 and t[2][2]==1.0:
            return True
        else:
            return False 
    def mapchains(self,c,newc):
        self.replicachainID_from_sourcechainID[c]=newc
        self.sourcechainID_from_replicachainID[newc]=c
    def get_replica_chainID(self,c):
        if c in self.replicachainID_from_sourcechainID:
           return self.replicachainID_from_sourcechainID[c]
        else:
            return c
    def report_chain_replicas(self):
        for k,v in self.replicachainID_from_sourcechainID.items():
            print('#### {} -> {}'.format(k,v))
    def get_base_chainID(self,newc):
        if newc in self.sourcechainID_from_replicachainID:
            return self.sourcechainID_from_replicachainID[newc]
        else:
            return newc
    def OneLiner(self):
        retstr=r'{ '
        for i in range(3):
            retstr+=r'{ '
            for j in range(4):
               retstr+='{} '.format(self.tmat[i][j])
            retstr+=r' } '
        retstr+='{ 0 0 0 1 } }'
        return retstr

