_PDBAtomNameDict_={'CL':'CLA'}
class Atom:
   def __init__(self,record_name,serial,name,altloc,resname,chainID,resseqnum,insertion,x,y,z,occ,beta,elem,charge):
      self.record_name=record_name.strip()
      self.serial=serial
      self.name=name.strip()
      self.altloc=altloc.strip()
      self.resname=resname.strip()
      self.chainID=chainID.strip()
      self.resseqnum=resseqnum
      self.insertion=insertion.strip()
      self.x=x
      self.y=y
      self.z=z
      self.occ=occ
      self.beta=beta
      self.elem=elem.strip()
      self.charge=charge.strip()
   def __str__(self):
      retstr='{}\n'+\
             '  serial    {:d}\n'+\
             '  name      {:s}\n'+\
             '  altloc    {:s}\n'+\
             '  resname   {:s}\n'+\
             '  chainID   {:s}\n'+\
             '  resseqnum {:d}\n'+\
             '  insertion {:s}\n'+\
             '  x         {:.3f}\n'+\
             '  y         {:.3f}\n'+\
             '  z         {:.3f}\n'+\
             '  occ       {:.2f}\n'+\
             '  beta      {:.2f}\n'+\
             '  elem      {:s}\n'+\
             '  charge    {:s}\n'
      return retstr.format(self.record_name,self.serial,self.name,self.altloc,self.resname,self.chainID,self.resseqnum,self.insertion,self.x,self.y,self.z,self.occ,self.beta,self.elem,self.charge)


