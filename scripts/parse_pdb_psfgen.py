''' Parses experimental PDB to 
       1. Extract loops of missing structure from REMARK 465 lines in a PDB
       2. Extract N-linked NAGs (more to come)
       3. Extract disulfides

    stdout can be redirected to a file that can be sourced by a mkpsf script
    after PDB is read in but before the 'guesscoords' statement

    Cameron F Abrams
    cfa22@drexel.edu
'''
_allowed_bond_types_=['DISU','NGLB','12xx','13xx','14xx','15xx','16xx','SA26E']
_allowed_ions_=['LIT','SOD','MG','POT','CAL','RUB','CES','BAR','ZN','CAD','CL']
_allowed_glycans_=['BMA','FUC','GAL','MAN','NAG','SIA']

_ResName123_={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY',
               'H':'HSE','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER',
               'T':'THR','W':'TRP','Y':'TYR','V':'VAL'}

_het_class_={'HOH':'WATER'}
_het_class_.update({k:'ION' for k in _allowed_ions_})
_het_class_.update({k:'GLYCAN' for k in _allowed_glycans_})
_het_class_.update({k:'PROTEIN' for k in _ResName123_.values()})
_ResDict_={'ZN':'ZN2','HOH':'TIP3','CL':'CLA','NAG':'BGNA','MAN':'AMAN','BMA':'BMAN','FUC':'AFUC','GAL':'BGAL','SIA':'ANE5AC'}
_AtNameDict_={'CL':'CLA'}

class bond:
   def __init__(self,ci,i,ai,cj,j,aj,t):
       self.ci=ci
       self.i=i
       self.ai=ai
       self.cj=cj
       self.j=j
       self.aj=aj
       self.type=t
       self.next=0
   def print_bonds(self):
       p=self
       while p!=0:
           if p.type in _allowed_bond_types_:
              cmd1=p.type[2]
              cmd2=p.type[3:]
              if p.type == 'NGLB':
                 p.cj=p.cj+'S'
                 print('patch {:2s}{}{} {}:{} {}:{}'.format(p.type[0:2],cmd1,cmd2,p.ci,p.i,p.cj,p.j))
              if p.type[0] == '1':
                 p.ci=p.ci+'S'
                 p.cj=p.cj+'S'
                 # ose_resid molid chain in_name
                 cmd1='[axeq {} 0 {} {} {}]'.format(p.i,p.ci[0],p.ai,p.j)
                 cmd2='[axeq {} 0 {} {} -1]'.format(p.j,p.cj[0],p.aj)
                 print('patch {:2s}{}{} {}:{} {}:{}'.format(p.type[0:2],cmd1,cmd2,p.cj,p.j,p.ci,p.i))
              if p.type == 'SA26E':
                 p.ci=p.ci+'S'
                 p.cj=p.cj+'S'
                 print('patch {:2s}{}{} {}:{} {}:{}'.format(p.type[0:2],cmd1,cmd2,p.ci,p.i,p.cj,p.j))
           else:
              print('ERROR: bond type {} not understood.'.format(p.type))
           p=p.next
   def add_bond(self,b):
       p=self
       while p.next!=0:
           p=p.next
       p.next=b

class mutopt:
   def __init__(self,lr,ri,rr,label):
       self.lr=lr
       self.ri=ri
       self.rr=rr
       self.label=label
       self.next=0
   def __str__(self):
       return '{}{}{} {}'.format(self.lr,self.ri,self.rr,self.label)
   def add_mutopt(self,nmut):
        p=self
        while p.next!=0:
           p=p.next
        p.next=nmut
   def print_inseg_option(self):
        print(r'   if {{ ${} == 1 }} {{'.format(self.label))
        print(r'       mutate {} {}'.format(self.ri,self.rr))
        print(r'   }')


_Missing_=0
_Present_=1
class loop:
   def __init__(self,state):
       self.chain=''
       self.rname=[]
       self.resid=[]
       self.next=0
       self.state=state
       self.mut=0
   def __str__(self):
       return 'chain {}, resid {} to {}, {}'.format(self.chain,self.resid[0],self.resid[-1],self.state)
   def add_res(self,c,r,i,s):
       if self.chain=='':
           self.index=0
           self.chain=c
           self.rname.append(r)
           self.resid.append(i)
           self.state=s
           return self
       else:
           if self.resid[-1]+1==i and self.state==s:
               self.rname.append(r)
               self.resid.append(i)
               return self
           elif i>self.resid[-1]+1 or self.chain!=c or self.state!=s:
#               print('add_res makes new loop {}'.format(c))
               newloop=loop(state=s)
               newloop.index=self.index+1
               newloop.chain=c
               newloop.rname.append(r)
               newloop.resid.append(i)
               self.next=newloop
               return newloop
           else:
               return self # unmolested
   def add_loop(self,me):
       p=self
       while p.next!=0:
          p=p.next
       p.next=me
   def remove(self,me):
       if self==me:
          self=self.next
          me.next=0
          return me,self
       else:
          p=self
          while p.next!=0:
             if p.next==me:
                p.next=me.next
                me.next=0
                return me,self
             p=p.next
#       print('cond 2')
       return 0,self
   def copy(self):
       retval=loop(self.state)
       retval.index=self.index
       retval.resid=self.resid[:]
       retval.rname=self.rname[:]
       retval.chain=self.chain
       retval.next=0
       return retval
   def pop_next_from_chain(self,c): # returns loop with lowest first resid
       min_resid=99999999
       p=self
       m=0
       while p!=0:
          if p.chain == c:
             if p.resid[0]<min_resid:
                m=p
                min_resid=p.resid[0]
          p=p.next
       if m==0:
#           print('cond 1')
           return 0,self
       return self.remove(m)
   def pop_first(self,state=_Missing_):
       p=self
       while p!=0:
          if p.state==state:
              return self.remove(p)
          p = p.next
       return 0,self
   def interleave(self):
       resorted=0
       while self!=0:
          l,self=self.pop()
          if resorted==0:
             resorted=l
             p=l
          else:
             p.next=l
             p=p.next
       self=resorted
       return self
   def make_chain_links(self):
       chains=set()
       p=self
       while p!=0:
           chains.add(p.chain)
           p=p.next
#       print(sorted(chains))
       chain_links=[[] for _ in range(len(chains))]
#       print(chain_links,len(chain_links),sorted(chains))
       for i,c in enumerate(sorted(chains)):
           p=self
           while p!=0:
               if p.chain==c:
                  x=p.copy()
                  if chain_links[i]==[]:
                     chain_links[i]=x
                  else:
                     chain_links[i].add_loop(x)
               p=p.next
       return chain_links
   def find_resid(self,r):
       p=self
       while p!=0:
           if r in p.resid:
               return p
           p = p.next
       return 0   
   def print_psfgen(self):
       p=self
       while p!=0:
           if p.state==_Present_:
              print(r'[atomselect top "chain {} and protein and resid {} to {}"] writepdb "{}_{}_to_{}.pdb"'.format(p.chain,p.resid[0],p.resid[-1],p.chain,p.resid[0],p.resid[-1]))
           p=p.next
       print('segment {} {{'.format(self.chain))
       p=self
       while p!=0:
           if p.state==_Present_:
               print('   pdb {}_{}_to_{}.pdb'.format(p.chain,p.resid[0],p.resid[-1]))
           elif p!=self and p.next!=0:
               for (r,i) in zip(p.rname,p.resid):
                   print('   residue {} {} {}'.format(i,r,p.chain))
           p = p.next
       p=self.mut
       while p!=0:
           p.print_inseg_option()
           p=p.next
       print('}')
       p=self
       while p!=0:
           if p.state==_Present_:
              print(r'coordpdb {}_{}_to_{}.pdb {}'.format(p.chain,p.resid[0],p.resid[-1],p.chain))
           p=p.next
       p=self
       while p.next!=0:
           if p.state==_Present_ and p.next.state==_Missing_ and p.next.next!=0:
              print(r'coord {} {} N [cacoIn_nOut {} {} 0]'.format(p.chain,p.next.resid[0],p.resid[-1],p.chain))
           p=p.next
   def print_loops(self,verbose=False):
        p=self
        while p != 0:
           if p.state==_Missing_:
              print('Missing Run #{}: chain {}, resid {} to {}'.format(p.index,p.chain,p.resid[0],p.resid[-1]))
           else:
              print('Present Run #{}: chain {}, resid {} to {}'.format(p.index,p.chain,p.resid[0],p.resid[-1]))
           if verbose:
               print('------------------------')
               for (r,i) in zip(p.rname,p.resid):
                   print('   residue {} {} {}'.format(i,r,p.chain))
               print('------------------------')
           p = p.next
   def get_chains(self):
       chains=set()
       p=self
       while p!=0:
          chains.add(p.chain)
          p=p.next
       return chains
   def print_loops_psfgen(self,state,ignore_free_ends=False):
       p=self
       lc=1
       while p!=0:
          if p.state==state:
             if ignore_free_ends==False or (p!=self and p.next!=0):
                print(r' {{ {} {:>4d} {:>4d} }} '.format(p.chain,p.resid[0],p.resid[-1]),end='')
                if lc%4==0:
                   print('\n           ',end='')
                lc+=1
          p=p.next

class glycres:
    def __init__(self,resname,resid):
       self.name=resname
       self.id=resid
       self.connections={}
    def add_connection(self,atom_name,dest_glycres,dest_atom_name):
       if not atom_name in self.connections:
          self.connections[atom_name]=(dest_glycres,dest_atom_name)
    def __str__(self,pad):
       s='{}{}-('.format(self.name,self.id)
       for ca,con in self.connections.items():
          s+='[{}]-[{}][{}]'.format(ca,con[1],str(con[0]))
       s+=')'
       return s

class glycan:
    def __init__(self,attach_resid):
       self.arid=attach_resid
       self.next=''
    def make_root(self,root_resid,root_resname):
       self.root=glycres(root_resname,root_resid)
    def add_res(self,at_res,at_atom,to_res,to_atom):
       at_res.add_connection(at_atom,to_res,to_atom)      

class het:
    def __init__(self,clas,nam,chn,rid):
       self.clas=clas
       self.chain=chn
       self.name=[nam]
       self.resid=[rid]
    def add_res(self,nam,rid):
       self.name.append(nam)
       self.resid.append(rid)
    def belongs_in_my_het(self,name,chain):
       return self.clas==_het_class_[name] and self.chain==chain

class hetlist:
    def __init__(self):
       self.h=[]
    # add a new HET either as a new segment/chain or as a residue in an existing chain
    def add_het(self,nam,chn,rid):
       q=''
       for h in self.h:
           if h.belongs_in_my_het(nam,chn):
              q=h
              break
       if q=='':
           self.h.append(het(_het_class_[nam],nam,chn,rid))
       else:
           q.add_res(nam,rid)
    def print_psfgen(self):
       for h in self.h:
#           print('HET name {} chain {}: resid {}-{}'.format(h.name,h.chain,min(h.resid),max(h.resid)))
           if h.clas == 'ION':
              segname=h.chain+'I'
           elif h.clas == 'GLYCAN':
              segname=h.chain+'S'
           elif h.clas == 'WATER':
              segname=h.chain+'WX'
           else:
              print('Error: unknown het class {}'.format(h.clas))
           segpdb='{}_{}_to_{}.pdb'.format(segname,min(h.resid),max(h.resid))
           print(r'set myseg [atomselect top "chain {} and resid {} to {}"]'.format(h.chain,min(h.resid),max(h.resid)))
           print(r'set sav_nm [$myseg get resname]')
           print(r'set new_nm [list]')   
           print(r'foreach r $sav_nm {')
           print(r'   lappend new_nm $RESDICT($r)')
           print(r'}')
           print(r'$myseg set resname $new_nm')
           print(r'set new_nm [list]')   
           print(r'set sav_nm [$myseg get name]')
           print(r'foreach r $sav_nm {')
           print(r'   if { [ info exists ANAMEDICT($r) ] } {')
           print(r'      lappend new_nm $ANAMEDICT($r)')
           print(r'   } else {')
           print(r'      lappend new_nm $r')
           print(r'   }')
           print(r'}')
           print(r'$myseg set name $new_nm')
           if h.clas=='WATER':
              print(r'$myseg set name OH2')
           print(r'$myseg writepdb "{}"'.format(segpdb))
           print(r'segment {} {{'.format(segname))
           print(r'   pdb {}'.format(segpdb))
           print(r'}')
           print(r'coordpdb {} {}'.format(segpdb,segname))

import sys


i=1
mut_csl=0
while i<len(sys.argv):
    if sys.argv[i]=='-pdb':
        i+=1
        fn=sys.argv[i]
    elif sys.argv[i]=='-mut':
        i+=1
        mut_csl=sys.argv[i]
    i+=1
M=0
if mut_csl!=0:
   muts=mut_csl.split(',')
   for m in muts:
      lr=_ResName123_[m[0]]
      rr=_ResName123_[m[-1]]
      ri=int(m[1:-1])
      if M==0:
         M=mutopt(lr,ri,rr,m)
      else:
         M.add_mutopt(mutopt(lr,ri,rr,m))
#m=M
#while m!=0:
#    print(m)
#    m=m.next

L=''
pl=L
H=hetlist()
D=''
#G=''
with open(fn) as f:
    for l in f:
        k=l.split()
        if k[0]=='COMPND' and k[1]=='3':
            chains=[]
            for i in range(3,len(k)):
                chains.append(k[i].strip(',').strip(';'))
        if k[0]=='REMARK' and k[1]=='465':
#            print(l)
            if len(k)==5:
               r=k[2]
               if r=='HIS':
                  r='HSE'
               c=k[3]
               i=int(k[4])
#               print(c,r,i)
               if L=='':
                   L=loop(_Missing_)
                   pl=L.add_res(c,r,i,_Missing_)
               else:
                   pl=pl.add_res(c,r,i,_Missing_)
        if k[0]=='ATOM':
            tr=l[17:20]
            if tr=='HIS':
               tr='HSE'
            tc=l[21:22]
            ti=int(l[22:26])
#            print(tc,ti)
            if L=='':
               L=loop(_Present_)
               pl=L.add_res(tc,tr,ti,_Present_)
            else:
#               print('adding {} {} {} present'.format(tc,tr,ti))
               pl=pl.add_res(tc,tr,ti,_Present_)

        if k[0]=='SSBOND':
            ci=l[15:16]
            i=int(l[16:22])
            cj=l[29:30]
            j=int(l[30:35])
            if D=='':
               D=bond(ci,i,'x',cj,j,'x','DISU')
            else:
               D.add_bond(bond(ci,i,'x',cj,j,'x','DISU'))
        if k[0]=='LINK':
            ai=l[10:15].strip()
            rni=l[17:20]
            ci=l[21:22]
            i=int(l[22:26])
            aj=l[40:45].strip()
            rnj=l[47:50]
            cj=l[51:52]
            j=int(l[52:56])
            ltyp='UNKNOWN'
            # for glycan bonds, C1 is always on the "i" resid
            if aj=='C1' and _het_class_[rni]=='GLYCAN':
               ai,aj=aj,ai
               rni,rnj=rnj,rni
               ci,cj=cj,ci
               i,j=j,i
            if rni=='ASN' and _het_class_[rnj]=='GLYCAN':
               ltyp='NGLB'
            elif _het_class_[rni]=='GLYCAN' and _het_class_[rnj]=='GLYCAN':
#               print(ai,aj)
               if ai=='C1' and aj=='O2':
                  ltyp='12xx'
               elif ai=='C1' and aj=='O3':
                  ltyp='13xx'
               elif ai=='C1' and aj=='O4':
                  ltyp='14xx'
               elif ai=='C1' and aj=='O5':
                  ltyp='15xx'
               elif ai=='C1' and aj=='O6':
                  ltyp='16xx'
#               elif ai=='C2' and aj=='O6':
#                  ltyp='SA26E'
               elif aj=='C2' and ai=='O6':
                  ltyp='SA26E'
            if D=='':
               D=bond(ci,i,ai,cj,j,aj,ltyp)
            else:
#               print('adding {}'.format(ltyp))
               D.add_bond(bond(ci,i,ai,cj,j,aj,ltyp))
        if k[0]=='HET':
            hetnam=l[7:10].strip()
            hetchain=l[12:13]
            hetresid=int(l[13:17])
            #print('DBG: hetnam {} hetchain {} hetresid {}'.format(hetnam,hetchain,hetresid))
            H.add_het(hetnam,hetchain,hetresid)
        if k[0]=='HETATM' and l[17:20]=='HOH':
            hetnam='HOH'
            hetchain=l[21:22]
            hetresid=int(l[22:26])
            H.add_het(hetnam,hetchain,hetresid)
#L.print_loops()
#print(L.get_chains())a

print(r'### autoPDBPARSE BEGIN')

for k,v in _ResDict_.items():
   print(r'set RESDICT({}) {}'.format(k,v))
for k,v in _AtNameDict_.items():
   print(r'set ANAMEDICT({}) {}'.format(k,v))

chains=sorted(L.get_chains())
cl=[0 for _ in range(len(chains))]
#print(chains)
i=0
for c in chains:
#   print('popping loops from chain',c)
   r,L=L.pop_next_from_chain(c)
#   print('popped',r,'remaining',L)
   if r!=0 and L==0:
       if cl[i]==0:
           cl[i]=r
       else:
           cl[i].add_loop(r)
   while r!=0 and L!=0:
      if cl[i]==0:
        cl[i]=r
      else:
        cl[i].add_loop(r)
#      print('popping loops from chain',c)
      r,L=L.pop_next_from_chain(c)
#      print('popped',r,'remaining',L)
   i+=1

#exit
#L.print_psfgen()
#G=D.build_glycans()
#print(cl)
print(r'set segs  {',end='')
for i,c in enumerate(cl):
   if i>0:
      print(r'           ',end='')
   c.print_loops_psfgen(_Present_)
   if i==len(cl)-1:
      print(' }')
   else:
      print()
print(r'set loops {',end='')
for i,c in enumerate(cl):
   if i>0:
      print(r'           ',end='')
   c.print_loops_psfgen(_Missing_,ignore_free_ends=True)
   if i==len(cl)-1:
      print(' }')
   else:
      print()
for c in cl:
   c.mut=M
   c.print_psfgen()
H.print_psfgen()
D.print_bonds()
