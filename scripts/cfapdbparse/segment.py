from residue import  _PDBResName123_, _pdb_glycans_, _pdb_ions_, _ResNameDict_PDB_to_CHARMM_, _ResNameDict_CHARMM_to_PDB_
_seg_class_={'HOH':'WATER'}
_seg_class_.update({k:'ION' for k in _pdb_ions_})
_seg_class_.update({k:'GLYCAN' for k in _pdb_glycans_})
_seg_class_.update({_ResNameDict_PDB_to_CHARMM_[k]:'GLYCAN' for k in _pdb_glycans_})
_seg_class_.update({k:'PROTEIN' for k in _PDBResName123_.values()})
_seg_class_['HIS']='PROTEIN'
_seg_class_['HSE']='PROTEIN'
_seg_class_['HSD']='PROTEIN'
_segname_second_character_={'PROTEIN':'','ION':'I','WATER':'W','GLYCAN':'G','OTHER':'O'}

class Run:
    ''' a set of contiguous residues with no gaps that can be loaded into a psfgen
        segment stanza by invoking a pdb file 
    '''
    def __init__(self,chainID,resseqnum1,resseqnum2):
        self.chainID=chainID
        self.resseqnum1=resseqnum1
        self.resseqnum2=resseqnum2
    def pdb_str(self):
        return '{}_{}_to_{}.pdb'.format(self.chainID,self.resseqnum1,self.resseqnum2)
    def __str__(self):
        return 'RUN: {} {} to {}'.format(self.chainID,self.resseqnum1,self.resseqnum2)

class Loop:
    ''' a set of contiguous residues from REMARK 465 pdb entries; i.e., they
        are missing from the base pdb file but present in the construct.  They
        must be included in a psfgen segment stanza via the 'residue' invocation.
    '''
    def __init__(self,chainID,resseqnum0,r):
        self.chainID=chainID
        self.resseqnum0=resseqnum0
        self.residues=[r]
        self.terminated=False
    def add_residue(self,r):
        self.residues.append(r)
    def __str__(self):
        return 'LOOP: {} ({})-[{} to {}]'.format(self.chainID,self.resseqnum0,self.residues[0].resseqnum,self.residues[-1].resseqnum)
    def caco_str(self):
        return 'coord {} {} N [cacoIn_nOut {} {} 0]\n'.format(self.chainID,self.residues[0].resseqnum,self.resseqnum0,self.chainID)

class Segment:
    def __init__(self,r,subcounter='',chain='',molid='top'):
        self.segname=r.chainID+_segname_second_character_[_seg_class_[r.name]]+subcounter
        self.source_chainID=r.source_chainID # survives cleavages
        self.parent_chain=chain
        self.segtype=_seg_class_[r.name]
        self.residues=[r]
        self.mutations=[]
        self.graft=''
        self.rootres=''
        self.molid=molid
        if _seg_class_[r.name]=='GLYCAN':
            #print(self.segname,r.name,r.resseqnum)
            self.rootres=r.up[0]
        r.segname=self.segname
        for a in r.atoms:
            a.segname=self.segname
    def __str__(self):
        retbase='({}){}[{}] {} - {}'.format(self.parent_chain.chainID,self.segname,self.segtype,self.residues[0].resseqnum,self.residues[-1].resseqnum)
        if self.rootres!='':
            r=self.rootres
            retbase='('+str(r)+')->'+retbase
        if self.graft!='':
            retbase+='::'+self.graft.source_pdb+':'+str(self.graft.source_segment)
        return retbase
    def add_residue(self,r):
        self.residues.append(r)
        r.segname=self.segname
        for a in r.atoms:
            a.segname=self.segname
    def apply_graft(self,g):
        self.graft=g
    def apply_attachment(self,a):
        self.attach=a
    def psfgen_segmentstanza(self):
        stanzastr='segment {} {{\n'.format(self.segname)
        if self.segtype=='PROTEIN':
            pdbs=[]
            Loops=[]
            Runs=[]
            cacostr=''
            makepdbstr=''
            inloop=False
            last_resseqnum=-1
            for r in self.residues:
                #print(r)
                if len(r.atoms)>0: # if this is not a missing resid
                    if last_resseqnum==-1: # initiate the list of runs with a new run
                        Runs.append(Run(r.chainID,r.resseqnum,-1))
                    elif inloop: # if we are in a loop
                        if len(Loops)>0: # terminate the latest loop
                            Loops[-1].terminated=True
                        inloop=False 
                        Runs.append(Run(r.chainID,r.resseqnum,-1)) # begin new run
                elif len(r.atoms)==0: # if this is a missing resid
                    if last_resseqnum==-1: # do not start segment with a loop of missing's
                        pass # Loops.append(Loop(r.chainID,-1,r))
                    elif not inloop: # terminate the current run and begin a new loop
                        Runs[-1].resseqnum2=last_resseqnum
                        Loops.append(Loop(r.chainID,last_resseqnum,r))
                    else: # we are in a loop, so add this to the latest loop
                        if len(Loops)>0:
                            Loops[-1].add_residue(r)
                    inloop=True
                last_resseqnum=r.resseqnum # ratchet up the tail resseqnum
            Runs[-1].resseqnum2=last_resseqnum
            #for r in Runs:
            #    print(r)
            #if len(Runs)==1: # if there is only one run
            # note, if we run out of residues while in a loop, the latest loop's 
            # 'terminated' flag is not set.  Later, we use this to prevent 
            # adding 'residue' statements to the psfgen segment stanza for
            # these residues
            # create stanza, mode 1: there ARE loops.  Generate pdb/residue statements for
            # each run/loop pair; but do nothing for the loop that is not terminated
            if len(Loops)>0:
                #print('### {} runs, {} loops'.format(len(Runs),len(Loops)))
                lim=max(len(Runs),len(Loops))
                for i in range(lim):
                    r='NULL'
                    l='NULL'
                    if i<len(Runs):
                       r=Runs[i]
                    if i<len(Loops):
                       l=Loops[i]
                    #print(r,l)
                    if r!='NULL':
                        thispdbstr,thispdb=self.psfgen_generate_input_pdb(r.resseqnum1,r.resseqnum2,r.pdb_str())
                        makepdbstr+=thispdbstr
                        pdbs.append(thispdb)
                        stanzastr+='   pdb {}\n'.format(thispdb)
                    if l!='NULL' and l.terminated==True:
                        for rr in l.residues:
                            if rr.name in _ResNameDict_PDB_to_CHARMM_:
                                nm=_ResNameDict_PDB_to_CHARMM_[rr.name]
                            else:
                                nm=rr.name
                            stanzastr+='   residue {}{} {} {}\n'.format(rr.resseqnum,rr.insertion,nm,rr.chainID)
            else: # this segment is a single run of residues that are all present in the pdb
                r=Runs[0]
                thispdbstr,thispdb=self.psfgen_generate_input_pdb(r.resseqnum1,r.resseqnum2,r.pdb_str())
                makepdbstr+=thispdbstr
                pdbs.append(thispdb)
                stanzastr+='   pdb {}\n'.format(thispdb)
            # specify mutations
            for m in self.mutations:
                stanzastr+=m.psfgen_segment_str()
            stanzastr+='}\n' # terminate the psfgen segment stanza
            loadcoordstr='# load/set coordinates from available data\n'
            for p in pdbs:
                loadcoordstr+='coordpdb {} {}\n'.format(p,self.segname)
            # for each loop, issue the CACO procecedure to set the coordinate of the N
            # in the first residue of the loop (guesscoord does not do this)
            if len(Loops)>0:
                r0=self.residues[0]
                for l in Loops:
                    if l.terminated:
                        cacostr+=l.caco_str()
            stanza=makepdbstr+stanzastr+loadcoordstr+cacostr
            return stanza,Loops

        elif self.segtype=='GLYCAN':
            r=Run(self.parent_chain.chainID,self.residues[0].resseqnum,self.residues[-1].resseqnum)
            pdb=r.pdb_str()
            makepdbstr,thispdb=self.psfgen_generate_input_pdb(r.resseqnum1,r.resseqnum2,pdb)
            stanzastr+='    pdb {}\n'.format(thispdb)
            stanzastr+='}\n'
            loadcoordstr='# load/set coordinates from available data\n'
            loadcoordstr+='coordpdb {} {}\n'.format(thispdb,self.segname)
            stanza=makepdbstr+stanzastr+loadcoordstr
            return stanza,[]
        else:
            r=Run(self.parent_chain.chainID,self.residues[0],self.residues[-1])
            pdb=r.pdb_str()
            makepdbstr,thispdb=self.psfgen_generate_input_pdb(r.resseqnum1,r.resseqnum2,pdb)
            stanzastr+='    pdb {}\n'.format(thispdb)
            stanzastr+='}\n'
            loadcoordstr='# load/set coordinates from available data\n'
            loadcoordstr+='coordpdb {} {}\n'.format(thispdb,self.segname)
            stanza=makepdbstr+stanzastr+loadcoordstr
            return stanza,[]
    def psfgen_generate_input_pdb(self,l,r,pdb_to_create):
        st=self.segtype
        molid=self.parent_chain.biomt.molid
        chainID=self.parent_chain.source_chainID
        g=self.graft
        a=self.attach
        p=pdb_to_create
        retstr=''
        if st=='PROTEIN':
            retstr+='set mysel [atomselect ${} "chain {} and resid {} to {}"]\n'.format(molid,chainID,l,r)
            retstr+='$mysel writepdb {}\n'.format(p)
        else:  # correct for any naming disagreements between PDB and CHARMM
            if g!='':
                g.ingraft_segname=self.segname
                g.ingraft_chainID=chainID
                retstr+='[atomselect ${} "chain {} and resid {} to {}"] writepdb {}\n'.format(molid,chainID,l,r,'GRAFTOVER'+p)
                retstr+='set ref [atomselect ${} "chain {} and resid {} and name C1 C2 O5 N2"]\n'.format(molid,g.target_chain,g.target_res)
                retstr+='set gra [atomselect {} "chain {} and resid {} and name C1 C2 O5 N2"]\n'.format(g.molid,g.source_chain,g.source_res1)
                retstr+='set refnum [$ref num]\n'
                retstr+='set granum [$gra num]\n'
                retstr+=r'if { $refnum != $granum } {'+'\n'
                retstr+='    $ref get resname\n'
                retstr+='    $gra get resname\n'
                retstr+='}\n'
                retstr+='set tm [measure fit $gra $ref]\n'
                retstr+='set myseg [atomselect {} "chain {} and resid {} to {}"]\n'.format(g.molid,g.source_chain,g.source_res1,g.source_res2)
                retstr+='set gra_orig_x [$myseg get x]\n'
                retstr+='set gra_orig_y [$myseg get y]\n'
                retstr+='set gra_orig_z [$myseg get z]\n'
                retstr+='$myseg move $tm\n'
                g.resid_dict={}
                for r in g.source_segment.residues:
                    g.resid_dict[r.resseqnum]=r.resseqnum+g.desired_offset
                retstr+='set savresid [$myseg get resid]\n'
                retstr+='set newresid [list]\n'
                retstr+=r'foreach oldresid [$myseg get resid] {'+'\n'
                retstr+='     lappend newresid [expr $oldresid + {:d}]\n'.format(g.desired_offset)
                retstr+='}\n'
                retstr+='$myseg set resid $newresid\n'
                #print('resid_dict:',g.resid_dict)
                p='{}_{}_to_{}-GRAFT.pdb'.format(g.source_chain,g.source_res1+g.desired_offset,g.source_res2+g.desired_offset)
            elif a!='':
                pass ### UNDER CONSTRUCTION -- ATTACH!!!
            else:
                retstr+='set myseg [atomselect ${} "chain {} and resid {} to {}"]\n'.format(molid,chainID,l,r)
            retstr+='set sav_resname [$myseg get resname]\n'
            retstr+='set new_resname [list]\n'
            retstr+=r'foreach r $sav_resname {'+'\n'
            retstr+=r'   if { [ info exists RESDICT($r) ] } {'+'\n'
            retstr+='      lappend new_resname $RESDICT($r)\n'
            retstr+=r'   } else {'+'\n'
            retstr+='      lappend new_resname $r\n'
            retstr+='   }\n'
            retstr+='}\n'
            retstr+='$myseg set resname $new_resname\n'
            retstr+='set new_name [list]\n'
            retstr+='set sav_name [$myseg get name]\n'
            retstr+=r'foreach r $sav_name {'+'\n'
            retstr+=r'   if { [ info exists ANAMEDICT($r) ] } {'+'\n'
            retstr+='      lappend new_name $ANAMEDICT($r)\n'
            retstr+=r'   } else {'+'\n'
            retstr+='      lappend new_name $r\n'
            retstr+='   }\n'
            retstr+='}\n'
            retstr+='$myseg set name $new_name\n'
            if st=='WATER':
                retstr+='$myseg set name OH2\n'
            
            retstr+='$myseg writepdb {}\n'.format(p)
            retstr+='# undo name changes and coordinate changes\n'
            retstr+='$myseg set resname $sav_resname\n'
            retstr+='$myseg set name $sav_name\n'
            if g!='' or a!='':
                retstr+='$myseg set resid $savresid\n'
                retstr+='$myseg set x $gra_orig_x\n'
                retstr+='$myseg set y $gra_orig_y\n'
                retstr+='$myseg set z $gra_orig_z\n'
        return retstr,p


