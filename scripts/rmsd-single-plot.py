import matplotlib.pyplot as plt
import argparse as ap
import numpy as np

def namd_config_to_dict(cfg):
    res={}
    with open(cfg,'r') as f:
        for l in f:
            tokens=l.lower().strip().split()
            #print(tokens)
            if len(tokens)>=2:
                res[tokens[0]]=tokens[1]
    return res

fac={}
p=ap.ArgumentParser()
p.add_argument('-rmsd',type=str,default=[],action='append',help='rmsd input file')
p.add_argument('-label',type=str,default=[],action='append',help='per-rmsd-file label')
p.add_argument('-namd-config',type=str,default=[],action='append',help='NAMD config file for info')
p.add_argument('-o',type=str,default='plot.png',help='output plot name')

args=p.parse_args()
for r in args.rmsd:
    fac[r]=1.0
xlabel='frame'
if len(args.label)==0:
    args.label=args.rmsd[:]
if len(args.namd_config)>0:
    namd_p={}
    for r,n in zip(args.rmsd,args.namd_config):
        namd_p[r]=namd_config_to_dict(n)
        fac[r]=float(namd_p[r]['timestep'])*float(namd_p[r]['dcdfreq'])/1.e6
    xlabel='time (ns)'

fig,ax=plt.subplots(1,1,figsize=(6,5))
ax.set_xlabel(xlabel)
ax.set_ylabel('RMSD (Å)')
t={}
f={}
for r,l in zip(args.rmsd,args.label):
    t[r],f[r]=np.loadtxt(r,unpack=True)
    ax.plot(t[r]*fac[r],f[r],label=l)
ax.legend()
plt.savefig(args.o,bbox_inches='tight')
