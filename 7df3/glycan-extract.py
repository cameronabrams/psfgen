#!/home/cfa/anaconda3/bin/python

with open('7df3.pdb','r') as fp:
   data=fp.read()
   lines=data.split('\n')


links=[]
for l in lines:
    if l[0:4]=='LINK':
        if l[13:16]=='ND2':
            links.append(l)

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

curr_maxresid_offset={}
print('[grafts]')
for l in links:
    n=int(l[22:26])
    gt=gtyp[n]
    gc=l[51:52]
    if gc not in curr_maxresid_offset:
        curr_maxresid_offset[gc]=1100
    else:
        curr_maxresid_offset[gc]+=100
    gr=int(l[52:56])
    outstr=ref[gt]+',{}:{},{}'.format(gc,gr,(0 if gr == 1 else gr+curr_maxresid_offset[gc]))
    print(outstr)
    

