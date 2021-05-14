import matplotlib.pyplot as plt
import argparse as ap
import numpy as np

def my_subst (target, vardict, vchar, delimchars='{}'):
    # entire target is a variable
    if target[0]==vchar and vchar+delimchars[0] not in target:
        return vardict[target[1:]]
    elif vchar in target:
        # there is a variable in there somewhere
        s=target.index(vchar)
        l=target.index(delimchars[0])
        r=target.index(delimchars[1])
        if l==s+1:
            repl=target[s:r+1]
            return target.replace(repl,vardict[target[l+1:r]])
    else:
        return target

def namd_config_to_dict(cfg):
    keys_w_multiple_values=['parameters']
    bad_chars = '; '
    res={}
    with open(cfg,'r') as f:
        for l in f:
            if l[0]=='#':
                continue
            tokens=l.lower().strip(bad_chars).split()
            for i in range(len(tokens)):
                tokens[i]=tokens[i].strip(bad_chars)
            if len(tokens)>=2:
                if tokens[0]=='set':
                    if 'set' not in res:
                        res['set']={}
                    res['set'][tokens[1]]=tokens[2]
                else:
                    if tokens[0] in keys_w_multiple_values:
                        if tokens[0] not in res:
                            res[tokens[0]]=[]
                        res[tokens[0]].append(tokens[1])
                    else:
                        res[tokens[0]]=tokens[1]
        for k,v in res.items():
            if k!='set':
                res[k]=my_subst(v,res['set'],'$')
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
        print('Temperature',namd_p[r]['langevintemp'])
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
