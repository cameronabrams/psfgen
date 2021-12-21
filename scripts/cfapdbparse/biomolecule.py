from moldata import MolData

""" biomolecule.py -- handling all biological assemblies in a PDB/mmCIF file """
class Biomolecule:
    ''' Container for handling info for "REMARK 350 BIOMOLECULE: #" stanzas in RCSB PDB files
        or _pdbx_struct blocks in mmCIF files '''
    def __init__(self,pdbrecord=None,cifdict=None,parent_molecule=''):
        self.chains_depot=[]
        self.biomt=[]
        self.pdbx_struct={}
        self.parent_molecule=parent_molecule
        self.apply_to_chainIDs=[]
        if pdbrecord==None and cifdict==None:
            # This constructor called with no pdbrecord and now cifdict 
            # creates the asymmetric unit, which is populated
            # by all explicit data in the input file
            self.index=0 # signifies asymmetric unit
            self.biomt.append(BiomT()) # give it the identity biomt
        elif pdbrecord!=None:
            if 'BIOMOLECULE:' in pdbrecord:
                # take the BIOMOLECULE index explicitly assigned in the PDB file
                try:
                    self.index=int(pdbrecord[23:25].strip())
                except:
                    print(f'Error parsing BIOMOLECULE statement: expected an integer and got {pdbrecord[23:25]}')
            else:
                print(f'Error parsing pdbrecord [{pdbrecord}]')
        elif cifdict!=None:
            self.index=int(cifdict['id'])
            for k in ['method_details','oligomeric_details','oligomeric_count']:
                self.pdbx_struct[k]=cifdict[k]
        else:
           print('Warning: Biomolecule.__init__() called with both a pdbrecord and cifdict; using pdbrecord')
            
    def parsePDBrecordwords(self,words):
        if len(words)>2:
            phrase=' '.join(words[2:])
            #print('biomolecule parse phrase',phrase)
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
                self.apply_to_chainIDs=[_.replace(',','') for _ in words[7:]]
            elif 'AND CHAINS:' in phrase:
                self.apply_to_chainIDs.extend([_.replace(',','') for _ in words[4:]])
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
            new_BiomT.apply_to_chainIDs=self.apply_to_chainIDs
            self.biomt.append(new_BiomT)
        self.biomt[-1].parseBIOMT(ax,words)

    def getBiomT(self,chainID=''):
        if len(chainID)>0:
            for t in self.biomt:
                if chainID in t.chainIDs:
                    return t
        return None
    
    def CIFBiomT(self,cifdict):
        self.biomt.append(BiomT(index=len(self.biomt)))
        self.biomt[-1].CIFBiomT(cifdict)

    ''' provide the moldata built for asymmetric unit during PDB readthru to the actual
        asymmetric unit'''
    def initializeAsymmetricUnit(self,md):
        b=self.biomt[0] # there is only one
        b.md=md
        apparent_chainIDs=list(set(b.md.Chains.keys()))
        for c in self.apply_to_chainIDs:
            if c not in apparent_chainIDs:
                print(f'PDB indicates AU has chain {c} but this chain is not represented among residues.')
        self.apply_to_chainIDs=apparent_chainIDs
        b.apply_to_chainIDs=self.apply_to_chainIDs
        #print(f'initializeAsymmetricUnit: {len(md.Mutations)} mutations')
    
    def inheritConstructs(self,au,cid_depot):
        for b in self.biomt:
            aumd=au.biomt[0].md
            ''' clone '''
#            for bcid in self.apply_to_chainIDs:
#               if bcid not in au.apply_to_chainIDs:
#                    print(f'PDB indicates Biomolecule {b.index} applies its transformations to chain {bcid} which is not in the AU')
            self.apply_to_chainIDs_as_read=self.apply_to_chainIDs[:]
            self.apply_to_chainIDs=au.apply_to_chainIDs[:]
            b.apply_to_chainIDs=self.apply_to_chainIDs
            b.apply_to_chainIDs_as_read=self.apply_to_chainIDs_as_read
            if b.isidentity():
                ''' this biomolecular assembly has one biomt that is the 
                    asymmetric unit '''
                b.md=au.biomt[0].md
            else:
                for aucid in au.apply_to_chainIDs:
                    try:
                        localcid=cid_depot.pop(0)
                    except:
                        print('Error: not enough alphabetical chain IDs!')
                        exit(-1)
                    b.mapChainIDs(aucid,localcid)
                
                b.md=aumd.Clone(chainmap=b.replicachainID_from_sourcechainID,invchainmap=b.sourcechainID_from_replicachainID)
                #print(f'inheritConstructs: {len(b.md.Mutations)} mutations.')
        return cid_depot                    

class BiomT:
    def __init__(self,index=0):
        self.md=None
        self.apply_to_chainIDs=[]
        self.apply_to_chainIDs_as_read=[]
#        self.ownChainIDs=[]
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
        print('{}BIOMT {:d} operates on chains {:s}'.format(indent,self.index,', '.join(self.apply_to_chainIDs)))
        if len(self.apply_to_chainIDs_as_read)>len(self.apply_to_chainIDs):
            print('{}{}(Note: This may not reflect the chains listed for this biomolecule in the PDB.'.format(indent))
            print('{}{}{}Original chains are {:s}'.format(indent,indent,indent,self.index,', '.join(self.apply_to_chainIDs_as_read)))
            print('{}{}One or more of these chains may have been deleted after processing mutations/deletions.)')
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

