from molecule import Molecule

def get_ra(mystr):
    if mystr[-1].isalpha():
       return int(mystr[:-1]),mystr[-1]
    else:
       return int(mystr),''

class Graft:
    def __init__(self,graftstr=''):
        # pdb,c:#-#,#,c:#
        self.graftstr=graftstr
        if len(graftstr)>0:
            dat=graftstr.split(',')
            if len(dat)==4:
                self.source_pdb=dat[0]
                source_chain_res=dat[1].split(':')
                self.source_rootres,self.source_rootins=get_ra(dat[2])
                target_chain_res=dat[3].split(':')
                if len(source_chain_res)==2:
                    self.source_chain=source_chain_res[0]
                    resrng=source_chain_res[1].split('-')
                    if len(resrng)==2:
                        for r in resrng:
                            self.source_res1,self.source_ins1=get_ra(r)
                    else:
                        print('ERROR: Malformed graft source resrange subsubargument: {}'.format(resrange))
                else:
                    print('ERROR: Malformed graft source subargument: {}'.format(source_chain_res))
                self.target_chain=target_chain_res[0]
                self.target_res,self.target_ins=get_ra(target_chain_res[1]) 
                self.molecule=Molecule(self.source_pdb)
                m=self.molecule
                for c in m.Chains:
                    c.MakeSegments(m.Links)
                    if c.chainID==self.source_chain:
                        for s in c.Segments:
                            for r in s.residues:
                                if r.resseqnum==self.source_res1:
                                    self.source_segment=s
                                    break
                if self.source_segment!='':
                    print('#### Graft {} will use this segment:'.format(self.graftstr))
                    print(self.source_segment)
                else:
                    print('ERROR: Could not find source segment for graft {}'.format(self.graftstr))         
            else:
               print('ERROR: Malformed graft argument: {}'.format(graftstr))
    def __str__(self):
        retstr='Graft from {:s}, chain {:s}, {:d}{} to {:d}{}, using {:d}{} as root, onto base chain {:s} {:d}{}'.format(self.source_pdb,self.source_chain,self.source_res1,self.source_ins1,self.source_res2,self.source_ins2,self.source_rootres,self.source_rootins,self.target_chain,self.target_res,self.target_ins)
        return retstr
    def load(self,fp,index):
        fp.write('mol new {}\n'.format(self.source_pdb))
        fp.write('set g{:d} [molinfo top get id]\n'.format(index))
        self.index=index

