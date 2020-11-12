class Biomolecule:
    ''' Container for handling info for "REMARK 350 BIOMOLECULE: #" stanzas in RCSB PDB files
        or _pdbx_struct blocks in RCSB CIF files '''
    def __init__(self,pdbrecord=None,cifdb=None):
        #print('__init__ with {}'.format(pdbrecord))
        self.chains=[]
        self.biomt=[]
        self.pdbx_struct={}
        if pdbrecord!=None and cifdb!=None:
           print('Warning: Biomolecule.__init__ called with both a pdbrecord and cifdb; using pdbrecord')
        if pdbrecord!=None:
            if 'BIOMOLECULE:' in pdbrecord:
               self.index=int(pdbrecord[23:25].strip())
           #print('new biomolecule index {:d}'.format(self.index))
        #self.author_biol_unit='None'
        #self.software_quat_struct='None'
        #self.software_used=['None']
        #self.total_buried_surface_area='None'
        #self.surface_area_units='None'
        #self.surface_area_complex='None'
        #self.change_solv_fe='None'
        #self.fe_units='None'
            if 'IDENTITY' in pdbrecord:  
                # caller would like an identity created, most likely because no 
                # BIOMOLECULE or PDBX_STRUCT was in the input file
                #print('#### identity biomt requested')
                self.index=1
                self.biomt.append(BiomT())
        elif cifdb!=None:
            self.index=int(cifdb['_pdbx_struct_assembly.id'])
            for k in ['method_details','oligomeric_details','oligomeric_count']:
                self.pdbx_struct[k]=cifdb['_pdbx_struct_assembly.{}'.format(k)]
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
    def show(self):
        print('Biomolecule {:d}'.format(self.index))
        print(self.pdbx_struct,self.chains)
        for b in self.biomt:
            b.show()
    def parseBIOMT(self,words):
        ax=int(words[2][-1])
        if ax==1:
            self.biomt.append(BiomT())
        self.biomt[-1].parseBIOMT(ax,words)

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
     def show(self):
         print('BIOMT {:d}'.format(self.index))
         print('    TMAT',self.tmat)
         print('    REPC',self.replicachainID_from_sourcechainID)
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

