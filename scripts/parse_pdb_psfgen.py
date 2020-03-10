''' Parses experimental PDB to 
       1. Extract loops of missing structure from REMARK 465 lines in a PDB
       2. Extract N-linked NAGs (more to come)
       3. Extract disulfides
    Cameron F Abrams
    cfa22@drexel.edu
'''
class ssbond:
   def __init__(self,ci,i,cj,j):
       self.ci=ci
       self.i=i
       self.cj=cj
       self.j=j
       self.next=0
   def print_ssbonds(self):
       p=self
       while p!=0:
           print('patch DISU {}:{} {}:{}'.format(p.ci,p.i,p.cj,p.j))
           p=p.next
   def add_ssbond(self,b):
       p=self
       while p.next!=0:
           p=p.next
       p.next=b
_Missing_=0
_Present_=1
class loop:
   def __init__(self,state):
       self.chain=''
       self.rname=[]
       self.resid=[]
       self.next=0
       self.state=state
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
       return 0,self
   def copy(self):
       retval=loop(self.state)
       retval.index=self.index
       retval.resid=self.resid[:]
       retval.rname=self.rname[:]
       retval.chain=self.chain
       retval.next=0
       return retval
   def pop_first(self,state=_Missing_):
       p=self
       while p!=0:
          if p.state==state:
              return self.remove(p)
          p = p.next
       return 0,self
   def interleave(self):
       ps=0
       fm=0
       fp=0
       while self!=0:
          fm,self=self.pop_first(state=_Missing_)
          if self!=0:
             fp,self=self.pop_first(state=_Present_)
          else:
             fp=0
          if ps==0:
             ps=fm
             ps.next=fp
             ps=ps.next
          else:
             ps.next=fm
             if fp!=0:
                ps=ps.next
                ps.next=fp
                ps=ps.next
       self=ps
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
              print('Missing Loop {}, chain {}: {} to {}'.format(p.index,p.chain,p.resid[0],p.resid[-1]))
           else:
              print('Present Segment {}, chain {}: {} to {}'.format(p.index,p.chain,p.resid[0],p.resid[-1]))
           if verbose:
               print('------------------------')
               for (r,i) in zip(p.rname,p.resid):
                   print('   residue {} {} {}'.format(i,r,p.chain))
               print('------------------------')

           p = p.next
   def print_loops_psfgen(self,state):
       p=self
       lc=1
       while p!=0:
          if p.state==state and p!=self and p.next!=0:
             print(r' {{ {} {:>4d} {:>4d} }} '.format(p.chain,p.resid[0],p.resid[-1]),end='')
             if lc%4==0:
                print('\n           ',end='')
             lc+=1
          p=p.next
import sys

fn=sys.argv[1]
L=loop(_Missing_)
pl=L
D=''
with open(fn) as f:
    for l in f:
        k=l.split()
        if k[0]=='COMPND' and k[1]=='3':
            chains=[]
            for i in range(3,len(k)):
                chains.append(k[i].strip(',').strip(';'))
        if k[0]=='REMARK' and k[1]=='465':
            if len(k)==5:
               r=k[2]
               if r=='HIS':
                  r='HSE'
               c=k[3]
               i=int(k[4])
#               print(c,r,i)
               pl=pl.add_res(c,r,i,_Missing_)
        if k[0]=='ATOM':
            tr=l[17:20]
            if tr=='HIS':
               tr='HSE'
            tc=l[21:22]
            ti=int(l[22:26])
#            print(tc,ti)
            pl=pl.add_res(tc,tr,ti,_Present_)
        if k[0]=='SSBOND':
            ci=l[15:16]
            i=int(l[16:22])
            cj=l[29:30]
            j=int(l[30:35])
            if D=='':
               D=ssbond(ci,i,cj,j)
            else:
               D.add_ssbond(ssbond(ci,i,cj,j))
#L.print_loops()
#L.print_psfgen()
cl=L.make_chain_links()
for c in cl:
   c.interleave()
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
   c.print_loops_psfgen(_Missing_)
   if i==len(cl)-1:
      print(' }')
   else:
      print()
for c in cl:
   c.print_psfgen()
#L.interleave()
D.print_ssbonds()
