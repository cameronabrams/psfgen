class BiomT_row:
    def __init__(self,ax,rep,a,b,c,d):
        self.ax=ax
        self.rep=rep-1 # make it start at zero
        self.row=[a,b,c,d]
    def __str__(self):
        retstr='BIOMT[{}] row {}: {} {} {} {}'.format(self.rep,self.ax,self.row[0],self.row[1],self.row[2],self.row[3])
        return retstr

class BiomT:
     def __init__(self,rows=[]):
         if len(rows)==3:
             self.rep=rows[0].rep
             self.tr=[rows[0].row,rows[1].row,rows[2].row]
         else:
             self.rep=0
             self.tr=[[1,0,0,0],[0,1,0,0],[0,0,1,0]]
         self.chainID_dict={}
         self.chainID_revdict={}
         self.molid='*'
         self.pdb='NONE'
     def __str__(self):
         retstr='BIOMT[{}]:\n'.format(self.rep)
         for i in range(len(self.tr)):
             for j in range(len(self.tr[i])):
                 retstr+='{:>10.4f} '.format(self.tr[i][j])
             retstr+='\n'
         return retstr
     def add_chain_replica(self,base_chainID,rep_chainID):
         self.chainID_dict[base_chainID]=rep_chainID
         self.chainID_revdict[rep_chainID]=base_chainID
     def get_replica_chainID(self,base_chainID):
         if base_chainID in self.chainID_dict:
             return self.chainID_dict[base_chainID]
         else:
             return '0'
     def get_base_chainID(self,rep_chainID):
         if rep_chainID in self.chainID_revdict:
             return self.chainID_revdict[rep_chainID]
         else:
             return '0'
     def OneLiner(self):
         retstr=r'{ '
         for i in range(3):
             retstr+=r'{ '
             for j in range(4):
                retstr+='{} '.format(self.tr[i][j])
             retstr+=r' } '
         retstr+='{ 0 0 0 1 } }'
         return retstr

