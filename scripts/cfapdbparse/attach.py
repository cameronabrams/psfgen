from molecule import Molecule

def get_ra(mystr):
    if mystr[-1].isalpha():
       return int(mystr[:-1]),mystr[-1]
    else:
       return int(mystr),''

def get_chain_ri1_ri2(mystr):
    s1=mystr.split(':')
    c=s1[0]
    s2=s1[1].split('-')
    r1,i1=get_ra(s2[0])
    if len(s2)==2:
        r2,i2=get_ra(s2[1])
    else:
        r2,i2='',''
    return c,r1,i1,r2,i2
    
class Attach:
    ''' Stores information to generate attachments.  An attachment refers to
        taking a segment or fragment from one pdb and inserting it in a specific
        way into the base molecule.  One must specify the source pdb, the chain
        and residue sequence of the atoms to be attached, which of those
        residues is the 'root', a chain and residue for the 'alignment basis', 
        and a target chain and residue to which the attachment is made.  
 
        The target residue and the alignment basis residue of the source must
        be congruent, since they serve to generate the transformation matrix
        to compute atomic coordinates.  (Currently, this condition is not checked.)
        The target residue and the root residue must constitute a pair for which
        a bond patch exists in the CHARMM force field.  Normally, the residue to which
        the piece to be transferred is attached is the same restype as the target residue,
        so it serves as the alignment basis.   Only the N, CA, C and O (backbone) serves
        as an alignment basis.
    '''
    def __init__(self,attachstr=''):
        # pdb,c:#-#,#,d:#,c:#,#
        self.attachstr=attachstr
        ''' pdb file containing the selection to be attached '''
        self.source_pdb=''
        ''' chain in that pdb file containing the selection to be attached '''
        self.source_chain=''
        ''' first residue of range of selection '''
        self.source_res1=''
        self.source_ins1=''
        ''' last residue of range of selection '''
        self.source_res2=''
        self.source_ins2=''
        ''' chain containing the residue used as an alignment basis '''
        self.source_align_chain=''
        ''' residue used as alignment basis '''
        self.source_align_res=''
        self.source_align_ins=''
        ''' segments in source molecule; built upon pdb read-in '''
        self.source_segment=''
        self.source_align_segment=''
        ''' chain in target (default is base molecule) to which attachment is made '''
        self.target_chain=''
        ''' target residue to which attachment is made '''
        self.target_res=''
        self.target_ins=''
        self.molecule=''
        self.molid=''
        self.index=''
        self.source_segment=''
        ''' resid offset for attached segment in target molecule '''
        self.desired_offset=''
        self.inattach_segname=''
        self.inattach_chainID=''
        self.resid_dict={}
        if len(attachstr)>0:
            dat=attachstr.split(',')
            if len(dat)==6:
                self.source_pdb=dat[0]
                self.source_chain,self.source_res1,self.source_ins1,self.source_res2,self.source_ins2=get_chain_ri1_ri2(dat[1])
                self.source_rootres,self.source_rootins=get_ra(dat[2])
                self.source_align_chain,self.source_align_res,self.source_align_ins,dum,dum=get_chain_ri1_ri2(dat[3])
                self.target_chain,self.target_res,self.target_ins,dum,dum=get_chain_ri1,di2(dat[4])
                self.desired_offset=int(dat[5])
                
                self.molecule=Molecule(self.source_pdb)
                m=self.molecule
                for c in m.Chains.values():
                    c.sort_residues()
                    c.MakeSegments(m.Links)
                for c in m.Chains.values():
                    if c.chainID==self.source_chain:
                        for s in c.Segments:
                            for r in s.residues:
                                if r.resseqnum==self.source_res1:
                                    self.source_segment=s
                                    break
                if self.source_segment!='':
                    print('#### Attachment {} will transfer this segment:'.format(self.graftstr))
                    print(self.source_segment)
                else:
                    print('ERROR: Could not find source segment for attachment {}'.format(self.graftstr))         
                for c in m.Chains.values():
                    if c.chainID==self.source_align_chain:
                        for s in c.Segments:
                            for r in s.residues:
                                if r.resseqnum==self.source_align_res:
                                    self.source_align_segment=s
                                    break
                if self.source_align_segment!='':
                    print('#### Attachment {} will align to target using residue {} of this segment:'.format(self.attachstr,self.source_align_res))
                    print(self.source_align_segment)
                else:
                    print('ERROR: Could not find source alignment segment for attachment {}'.format(self.attachstr))         

            else:
               print('ERROR: Malformed attach argument: {}'.format(attachstr))
    def attachStr(self,replace_targ_chain=''):
        ''' regenerate the string code for this attachment with option to change target chain ID '''
        # pdb,c:#-#,#,d:#,c:#,#
        return '{},{}:{}{}-{}{},{}{},{}:{}{},{}:{}{},{}'.format(self.source_pdb,self.source_chain,self.source_res1,self.source_ins1,self.source_res2,self.source_ins2,self.source_rootres,self.source_rootins,replace_targ_chain if replace_targ_chain != '' else self.target_chain,self.target_res,self.target_ins,self.desired_offset)
    def __str__(self):
        retstr='Attach from {:s}, chain {:s}, {:d}{} to {:d}{}, using {:d}{} as root and chain {:s} resid {:d}{} as an alignment basis, onto base chain {:s} {:d}{} with resid offset {:d}'.format(self.source_pdb,self.source_chain,self.source_res1,self.source_ins1,self.source_res2,self.source_ins2,self.source_rootres,self.source_rootins,self.target_chain,self.target_res,self.target_ins,self.desired_offset)
        return retstr

    def load(self,fp,index):
        fp.write('mol new {}\n'.format(self.source_pdb))
        fp.write('set g{:d} [molinfo top get id]\n'.format(index))
        self.molid='$g{:d}'.format(index)
        self.index=index
   
