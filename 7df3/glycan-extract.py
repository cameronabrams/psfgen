#!/home/cfa/anaconda3/bin/python
# 
# This script generates a modfile for cfapdbparse.py
# that instructs it to graft glycans onto the NAGs
# Each glycan is one of 'complex', 'hybrid', or
# 'oligomannose'.  These are set explicitly below
# in the `gtyp` dict.  The `ref` dict indicates
# the source pdbs for each type; each value in this
# dictionary is of the form
#
# [PDBfile],[Chain]:[Resid-range],[Base-resid]
#
# where [PDBfile] is the source PDB file for the glycan,
# [Chain] is the chain of the glycan, [Resid-range]
# is the range of resids for the glycan, and [Base-resid]
# is the residue that is used as an alignment basis when
# grafting the glycan onto an existing glycan stub (typically
# one or two NAGs in the RCSB file).
#
# Cameron F Abrams cfa22@drexel.edu
#
import argparse as ap
parser=ap.ArgumentParser()
parser.add_argument('-pdb',type=str,default='',help='RCSB pdb file')
parser.add_argument('-o','--outfile',type=str,default='glycans.mod',help='modsfile to generate')
args=parser.parse_args()

# borrowed from Watanabe 2020 fof SARS-CoV-2 S
gtyp={}
ref={}
gtyp[17]='c'
gtyp[61]='h'
gtyp[122]='h'
gtyp[149]='c'
gtyp[165]='c'
gtyp[234]='o'
gtyp[282]='c'
gtyp[331]='c'
gtyp[343]='c'
gtyp[603]='h'
gtyp[616]='c'
gtyp[657]='c'
gtyp[709]='o'
gtyp[717]='h'
gtyp[801]='h'
gtyp[1074]='h'
gtyp[1098]='c'
gtyp[1134]='c'
ref['o']='2wah.pdb,C:1-9,1'
ref['h']='4b7i.pdb,C:1-8,1'
ref['c']='4byh.pdb,C:1-10,1'

# read entire pdb file
with open(args.pdb,'r') as fp:
   data=fp.read()
   lines=data.split('\n')

# establish the list of LINKs that connect to ND2 atoms
# -- these are N-linked glycans
links=[]
for l in lines:
    if l[0:4]=='LINK':
        if l[13:16]=='ND2':
            links.append(l)

# parse each link into a modsfileline in a [grafts] section
fp=open(args.outfile,'w')
curr_maxresid_offset={}
fp.write('[grafts]\n')
for l in links:
    n=int(l[22:26])
    gt=gtyp[n]
    gc=l[51:52]
    if gc not in curr_maxresid_offset:
        curr_maxresid_offset[gc]=1100
    else:
        curr_maxresid_offset[gc]+=100
    gr=int(l[52:56])
    outstr=ref[gt]+',{}:{},{}\n'.format(gc,gr,(0 if gr == 1 else gr+curr_maxresid_offset[gc]))
    fp.write(outstr)
fp.close()
print('Generated {} from LINKs in {}.'.format(args.outfile,args.pdb))
