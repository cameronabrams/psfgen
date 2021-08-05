'''functions for generating TcL/VMD atomselect commands'''

def backup(selname):
    retstr=''
    attr=['x','y','z','resid','resname','name']
    for a in attr:
       retstr+='set {}_orig_{} [${} get {}]\n'.format(selname,a,selname,a)
    return retstr

def restore(selname):
    retstr=''
    attr=['x','y','z','resid','resname','name']
    for a in attr:
       retstr+='${} set {} ${}_orig_{}\n'.format(selname,a,selname,a)
    return retstr

def residshift(selname,shift):
    retstr=''
    retstr+='set new_resid [list]\n'
    retstr+=r'foreach oldresid [${} get resid] {{'.format(selname)+'\n'
    retstr+='     lappend new_resid [expr $oldresid + {:d}]\n'.format(shift)
    retstr+='}\n'
    retstr+='${} set resid $new_resid\n'.format(selname)
    return retstr

''' charmm_namify converts commonly found atom and residue names in PDB files to their 
    appropriate charmm names -- this is only used for NON-PROTEIN SEGMENTS.  In atom.py, you
    must put an associated entry in the global _PDBAtomNameDict_ in order for this to work. '''
def charmm_namify(selname,iswater=False):
    retstr=''
    retstr+='set new_resname [list]\n'
    retstr+=r'foreach r [${} get resname] {{'.format(selname)+'\n'
    retstr+=r'   if { [ info exists RESDICT($r) ] } {'+'\n'
    retstr+='      lappend new_resname $RESDICT($r)\n'
    retstr+=r'   } else {'+'\n'
    retstr+='      lappend new_resname $r\n'
    retstr+='   }\n'
    retstr+='}\n'
    retstr+='set new_name [list]\n'
    retstr+=r'foreach r [${} get name] {{'.format(selname)+'\n'
    retstr+=r'   if { [ info exists ANAMEDICT($r) ] } {'+'\n'
    retstr+='      lappend new_name $ANAMEDICT($r)\n'
    retstr+=r'   } else {'+'\n'
    retstr+='      lappend new_name $r\n'
    retstr+='   }\n'
    retstr+='}\n'
    retstr+='${} set resname $new_resname\n'.format(selname)
    retstr+='${} set name $new_name\n'.format(selname)
    if iswater=='WATER':
        retstr+='${} set name OH2\n'.format(selname)
    return retstr
