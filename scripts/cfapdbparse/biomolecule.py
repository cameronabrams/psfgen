""" biomolecule.py -- handling all biological assemblies in a PDB/mmCIF file """
class Biomolecule:
    ''' Container for handling info for "REMARK 350 BIOMOLECULE: #" stanzas in RCSB PDB files
        or _pdbx_struct blocks in mmCIF files '''
    def __init__(self,pdbrecord=None,cifdict=None):
        self.chains_depot=[]
        self.biomt=[]
        self.pdbx_struct={}
        if pdbrecord==None and cifdict==None:
            # this is interpreted as the "asymmetric unit"; i.e., an assembly composed of only those
            # atoms explicitly in the PDB/mmCIF file
            self.index=0 # signifies asymmetric unit
            self.biomt.append(BiomT()) # give it the identity biomt
        elif pdbrecord!=None:
            if 'BIOMOLECULE:' in pdbrecord:
                self.index=int(pdbrecord[23:25].strip())  # take the index explicitly assigned in the PDB file
            else:
                print('Error:  Cannot parse PDBRecord [{}]'.format(pdbrecord))
        elif cifdict!=None:
            self.index=int(cifdict['id'])
            for k in ['method_details','oligomeric_details','oligomeric_count']:
                self.pdbx_struct[k]=cifdict[k]
        else:
           print('Warning: Biomolecule.__init__() called with both a pdbrecord and cifdict; using pdbrecord')

            
    def parsePDBrecordwords(self,words):
        if len(words)>2:
            phrase=' '.join(words[2:])
          #  print(phrase)
            if 'AUTHOR DETERMINED BIOLOGICAL UNIT:' in phrase:
           #     print(words[-1])
                self.pdbx_struct['author_biol_unit']=words[-1]
#                self.pdbx_struct['software_used']=[_.replace(',','') for _ in words[4:]]
            elif 'SOFTWARE DETERMINED QUATERNARY STRUCTURE:' in phrase:
           #     print(words[-1])
                self.pdbx_struct['software_biol_unit']=words[-1]
#                self.pdbx_struct['software_used']=[_.replace(',','') for _ in words[4:]]
            elif 'SOFTWARE USED:' in phrase:
           #     print(words[-1])
                self.pdbx_struct['software_used']=words[-1]
#                self.pdbx_struct['software_used']=[_.replace(',','') for _ in words[4:]]
            elif 'TOTAL BURIED SURFACE AREA:' in phrase:
                self.pdbx_struct['total_buried_surface_area']=int(words[-2])
                self.pdbx_struct['surface_area_units']=words[-1]
            elif 'SURFACE AREA OF THE COMPLEX' in phrase:
                self.pdbx_struct['surface_area_complex']=int(words[-2])
            elif 'CHANGE IN SOLVENT FREE ENERGY:' in phrase:
                self.pdbx_struct['change_solv_fe']=float(words[-2])
                self.pdbx_struct['fe_units']=words[-1]
            elif 'APPLY THE FOLLOWING TO CHAINS:' in phrase:
                self.chain_depot=[_.replace(',','') for _ in words[7:]]
            elif 'AND CHAINS:' in phrase:
                self.chain_depot.extend([_.replace(',','') for _ in words[4:]])
            elif 'BIOMT' in words[2]:
                self.parseBIOMT(words)
            else:
                print('#### ERROR: Unrecognized PDB-format REMARK 350 line:')
                print(' '.join(words))
    
    def show(self,isActive=False,indent='    '):
        if self.index==0:
            desig='Asymmetric Unit'
        else:
            desig='Biological Assembly'
        activeLabel='**ACTIVE**' if isActive else ''
        print('{}{} {} {:d} {}'.format(indent,activeLabel,desig,self.index,activeLabel))
        if len(self.pdbx_struct)>0:
            print(indent*2,self.pdbx_struct)
        for b in self.biomt:
            b.show(indent*2)
    
    def parseBIOMT(self,words):
        ax=int(words[2][-1])
        if ax==1:   # "BIOMT1" detected
            new_BiomT=BiomT(index=len(self.biomt))
            # empty the chain_depot
            if len(self.chain_depot)==0:
                print('Warning: New BIOMT detected but no chains are designated for it.')
            else:
                new_BiomT.chainIDs=self.chain_depot[:]
                self.chain_depot=[]
            self.biomt.append(new_BiomT)
        self.biomt[-1].parseBIOMT(ax,words)
    
    def CIFBiomT(self,cifdict):
        self.biomt.append(BiomT(index=len(self.biomt)))
        self.biomt[-1].CIFBiomT(cifdict)

class BiomT:
    def __init__(self,index=0):
        self.chainIDs=[]
        self.index=index
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
        print('{}BIOMT {:d} operates on chains {:s}'.format(indent,self.index,', '.join(self.chainIDs)))
        if not self.isidentity():
            print('{}    TMAT'.format(indent),self.tmat)
            if len(self.replicachainID_from_sourcechainID)>0:
                print('{}    REPC'.format(indent),self.replicachainID_from_sourcechainID)
        else:
            print('{}    IDENTITY'.format(indent))
    def isidentity(self):
        t=self.tmat
        if t[0][0]==1.0 and t[1][1]==1.0 and t[2][2]==1.0:
            return True
        else:
            return False 
    def mapChainIDs(self,c,newc):
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

