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
    def __init__(self,chainID,resseqnum1,resseqnum2):
        self.chainID=chainID
        self.resseqnum1=resseqnum1
        self.resseqnum2=resseqnum2
    def pdb_str(self):
        return '{}_{}_to_{}.pdb'.format(self.chainID,self.resseqnum1,self.resseqnum2)
    def __str__(self):
        return 'RUN: {} {} to {}'.format(self.chainID,self.resseqnum1,self.resseqnum2)

class Loop:
    def __init__(self,chainID,resseqnum0,r):
        self.chainID=chainID
        self.resseqnum0=resseqnum0
        self.residues=[r]
        self.terminated=False
    def add_residue(self,r):
        self.residues.append(r)
    def __str__(self):
        return 'LOOP chainID {} r0 {} r1 {} to r2 {}'.format(self.chainID,self.resseqnum0,self.residues[0].resseqnum,self.residues[-1].resseqnum)
    def caco_str(self):
        return 'coord {} {} N [cacoIn_nOut {} {} 0]\n'.format(self.chainID,self.residues[0].resseqnum,self.resseqnum0,self.chainID)

class Segment:
    def __init__(self,r,subcounter=''):
        self.segname=r.chainID+_segname_second_character_[_seg_class_[r.name]]+subcounter
        self.source_chainID=r.source_chainID # survives cleavages
        self.segtype=_seg_class_[r.name]
        self.residues=[r]
        self.mutations=[]
        r.segname=self.segname
        for a in r.atoms:
            a.segname=self.segname
    def __str__(self):
        return '{}[{}] {} - {}'.format(self.segname,self.segtype,self.residues[0].resseqnum,self.residues[-1].resseqnum)
    def add_residue(self,r):
        self.residues.append(r)
        r.segname=self.segname
        for a in r.atoms:
            a.segname=self.segname
    def psfgen_segmentstanza(self):
        pdbs=[]
        Loops=[]
        Runs=[]
        cacostr=''
        suppstr=''
        retstr='segment {} {{\n'.format(self.segname)
        inloop=False
        last_resseqnum=-1
#        for r in self.residues:
#            fp.write(r)
        for r in self.residues:
            if len(r.atoms)>0:
                if last_resseqnum==-1:
#                    fp.write('appending',r)
                    Runs.append(Run(r.chainID,r.resseqnum,-1))
                elif inloop:
                    if len(Loops)>0:
                        Loops[-1].terminated=True
                    inloop=False
                    Runs.append(Run(r.chainID,r.resseqnum,-1))
                else:
#                    fp.write('2passing',r)
                    pass
            elif len(r.atoms)==0:
                if last_resseqnum==-1:
                    pass # Loops.append(Loop(r.chainID,-1,r))                    
                elif not inloop:
                    Runs[-1].resseqnum2=last_resseqnum
                    Loops.append(Loop(r.chainID,last_resseqnum,r))
                else:
#                    fp.write(last_resseqnum)
                    if len(Loops)>0:
                        Loops[-1].add_residue(r)
                    else:
                        pass
                inloop=True
            else:
#                fp.write('passing', r)
                pass
            last_resseqnum=r.resseqnum
#        fp.write(self.segname,len(Runs))
        if len(Runs)==1:
            Runs[0].resseqnum2=last_resseqnum
#            fp.write('END OF LOOP: segname {} r.resseqnum {} {} last_resseqnum {}'.format(self.segname,r.resseqnum,r.name,last_resseqnum))
        if len(Loops)>0:
            for r,l in zip(Runs,Loops):
                pdbs.append(r.pdb_str())
                suppstr+=psfgen_write_pdb(self.segtype,self.source_chainID,r.chainID,r.resseqnum1,r.resseqnum2,pdbs[-1])
                retstr+='   pdb {}\n'.format(pdbs[-1])
                if l.terminated==True:
                    for rr in l.residues:
                        if rr.name in _ResNameDict_PDB_to_CHARMM_:
                            nm=_ResNameDict_PDB_to_CHARMM_[rr.name]
                        else:
                            nm=rr.name
                        retstr+='   residue {} {} {}\n'.format(rr.resseqnum,nm,rr.chainID)
            for m in self.mutations:
                retstr+=m.psfgen_segment_str()
        else:
            r=Runs[0]
            pdbs.append(r.pdb_str())
            suppstr+=psfgen_write_pdb(self.segtype,self.source_chainID,r.chainID,r.resseqnum1,r.resseqnum2,pdbs[-1])
            retstr+='   pdb {}\n'.format(pdbs[-1])

        retstr+='}'
        coordstr=''
        for p in pdbs:
            coordstr+='coordpdb {} {}\n'.format(p,self.segname)
        if len(Loops)>0:
            r0=self.residues[0]
            for l in Loops:
                if l.terminated:
                    cacostr+=l.caco_str()
        return retstr,suppstr,coordstr,cacostr,Loops

def psfgen_write_pdb(st,source,c,l,r,p):
    if st=='PROTEIN':
        retstr='[atomselect top "chain {} and resid {} to {}"] writepdb {}\n'.format(source,l,r,p)
    else:
        retstr='# Changing resnames and atom names in {}\n'.format(p)
        retstr+='set myseg [atomselect top "chain {} and resid {} to {}"]\n'.format(source,l,r,p)
        retstr+='set sav_nm [$myseg get resname]\n'
        retstr+='set new_nm [list]\n'
        retstr+=r'foreach r $sav_nm {'+'\n'
        retstr+=r'   if { [ info exists RESDICT($r) ] } {'+'\n'
        retstr+='      lappend new_nm $RESDICT($r)\n'
        retstr+=r'   } else {'+'\n'
        retstr+='      lappend new_nm $r\n'
        retstr+='   }\n'
        retstr+='}\n'
        retstr+='$myseg set resname $new_nm\n'
        retstr+='set new_nm [list]\n'
        retstr+='set sav_nm [$myseg get name]\n'
        retstr+=r'foreach r $sav_nm {'+'\n'
        retstr+=r'   if { [ info exists ANAMEDICT($r) ] } {'+'\n'
        retstr+='      lappend new_nm $ANAMEDICT($r)\n'
        retstr+=r'   } else {'+'\n'
        retstr+='      lappend new_nm $r\n'
        retstr+='   }\n'
        retstr+='}\n'
        retstr+='$myseg set name $new_nm\n'
        if st=='WATER':
            retstr+='$myseg set name OH2\n'
        retstr+='$myseg writepdb {}\n'.format(p)

    return retstr


