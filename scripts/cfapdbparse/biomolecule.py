class Biomolecule:
    ''' Container for handling info for "REMARK 350 BIOMOLECULE: #" stanzas in RCSB PDB files
    '''
    def __init__(self,pdbrecord):
       tok=pdbrecord[11:23]
       #print('__init__ with {}'.format(pdbrecord))
       if tok=='BIOMOLECULE:':
           self.index=int(pdbrecord[23:25].strip())
           #print('new biomolecule index {:d}'.format(self.index))
           self.chains=[]
           self.biomt=[]
           self.author_biol_unit='None'
           self.software_quat_struct='None'
           self.software_used=['None']
           self.total_buried_surface_area='None'
           self.surface_area_units='None'
           self.surface_area_complex='None'
           self.change_solv_fe='None'
           self.fe_units='None'
           self.chains=[]
    def parsePDBrecordwords(self,words):
        if len(words)>2:
            phrase=' '.join(words[2:])
          #  print(phrase)
            if 'AUTHOR DETERMINED BIOLOGICAL UNIT:' in phrase:
           #     print(words[-1])
                self.author_biol_unit=words[-1]
            elif 'SOFTWARE DETERMINED QUATERNARY STRUCTURE:' in phrase:
                self.software_quat_struct=words[-1]
            elif 'SOFTWARE USED:' in phrase:
                self.software_used=[_.replace(',','') for _ in words[4:]]
            elif 'TOTAL BURIED SURFACE AREA:' in phrase:
                self.total_buried_surface_area=int(words[-2])
                self.surface_area_units=words[-1]
            elif 'SURFACE AREA OF THE COMPLEX' in phrase:
                self.surface_area_complex=int(words[-2])
            elif 'CHANGE IN SOLVENT FREE ENERGY:' in phrase:
                self.change_solv_fe=float(words[-2])
                self.fe_units=words[-1]
            elif 'APPLY THE FOLLOWING TO CHAINS:' in phrase:
                self.chains=[_.replace(',','') for _ in words[7:]]
            elif 'BIOMT' in words[2]:
                self.parseBIOMT(words)
            else:
                print('#### ERROR: Unrecognized REMARK 350 line:')
                print(' '.join(words))
    def show(self):
        print('Biomolecule {:d}'.format(self.index))
        print(self.author_biol_unit,self.software_quat_struct,self.total_buried_surface_area,self.surface_area_complex,self.change_solv_fe,self.chains)
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
         self.tmat=[[],[],[]]
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
         print(self.tmat)
         print(self.replicachainID_from_sourcechainID)
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

