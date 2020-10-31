from molecule import Molecule
import sel

def get_ra(mystr):
    if mystr[-1].isalpha():
       return int(mystr[:-1]),mystr[-1]
    else:
       return int(mystr),''

class Graft:
    def __init__(self,graftstr=''):
        # pdb,c:#-#,#,c:#,#
        self.graftstr=graftstr
        self.source_pdb=''
        self.source_chain=''
        self.source_res1=''
        self.source_ins1=''    
        self.source_res2=''
        self.source_ins2=''
        self.source_rootres=''
        self.source_rootins=''
        self.target_chain=''
        self.target_res=''
        self.target_ins=''
        self.source_segment=''
        self.source_segment=''
        self.desired_offset=''
        self.transformed_pdb=''
        if len(graftstr)>0:
            dat=graftstr.split(',')
            if len(dat)==5:
                self.source_pdb=dat[0]
                source_chain_res=dat[1].split(':')
                self.source_rootres,self.source_rootins=get_ra(dat[2])
                target_chain_res=dat[3].split(':')
                self.desired_offset=int(dat[4])
                if len(source_chain_res)==2:
                    self.source_chain=source_chain_res[0]
                    resrng=source_chain_res[1].split('-')
                    if len(resrng)==2:
                        self.source_res1,self.source_ins1=get_ra(resrng[0])
                        self.source_res2,self.source_ins2=get_ra(resrng[1])
                    else:
                        print('ERROR: Malformed graft source resrange subsubargument: {}'.format(resrange))
                else:
                    print('ERROR: Malformed graft source subargument: {}'.format(source_chain_res))
                self.target_chain=target_chain_res[0]
                self.target_res,self.target_ins=get_ra(target_chain_res[1]) 
                self.molecule=Molecule(self.source_pdb,isgraft=True)
                m=self.molecule
                for c in m.Chains.values():
                    c.sort_residues()
                    c.MakeSegments(m.Links)
                if self.source_chain in m.Chains:
                    c=m.Chains[self.source_chain]
                    for s in c.Segments:
                        for r in s.residues:
                            if r.resseqnum==self.source_res1:
                               self.source_segment=s
                               break
                #if self.source_segment!='':
                if self.source_segment=='':
                #    print('#### Graft {} will copy from segment {}'.format(self.graftstr.strip(),str(self.source_segment)))
                #else:
                    print('ERROR: Could not find source segment for graft {}'.format(self.graftstr))         
            else:
                print('ERROR: Malformed graft argument: {}'.format(graftstr))
    def graftStr(self,replace_targ_chain=''):
        return '{},{}:{}{}-{}{},{}{},{}:{}{},{}'.format(self.source_pdb,self.source_chain,self.source_res1,self.source_ins1,self.source_res2,self.source_ins2,self.source_rootres,self.source_rootins,replace_targ_chain if replace_targ_chain != '' else self.target_chain,self.target_res,self.target_ins,self.desired_offset)
    def __str__(self):
        retstr='Graft from {:s}, chain {:s}, {:d}{} to {:d}{}, using {:d}{} as root, onto base chain {:s} {:d}{} with resid offset {:d}'.format(self.source_pdb,self.source_chain,self.source_res1,self.source_ins1,self.source_res2,self.source_ins2,self.source_rootres,self.source_rootins,self.target_chain,self.target_res,self.target_ins,self.desired_offset)
        return retstr

    def load(self,fp):
        self.molecule.load(fp)

    def transform(self,base):
        a=''
        a+='set ref [atomselect ${} "chain {} and resid {} and noh"]\n'.format(base.molid_varname,self.target_chain,self.target_res)
        a+='set gra [atomselect ${} "chain {} and resid {} and noh"]\n'.format(self.molecule.molid_varname,self.source_chain,self.source_res1)
        a+='set tra [atomselect ${} "chain {} and resid {} to {}"]\n'.format(self.molecule.molid_varname,self.source_chain,self.source_res1,self.source_res2)
        a+=r'if { [$ref num] != [$gra num] } {'+'\n'
        a+='    puts "psfgen: warning: target and graft alignment references are not congruent"\n'
        a+='}\n'
        a+=sel.backup('tra')
        a+='$tra move [measure fit $gra $ref]\n'
        self.resid_dict={}
        for r in self.source_segment.residues:
            self.resid_dict[r.resseqnum]=r.resseqnum+self.desired_offset
        a+=sel.residshift('tra',self.desired_offset)
        self.transformed_pdb='{}-{}_{}_to_{}-GRAFT.pdb'.format(self.ingraft_chainID,self.source_chain,self.source_res1+self.desired_offset,self.source_res2+self.desired_offset)
        a+=sel.charmm_namify('tra')
        a+='$tra writepdb {}\n'.format(self.transformed_pdb)
        a+=sel.restore('tra')
        self.transform_commands=a
        return a
 
