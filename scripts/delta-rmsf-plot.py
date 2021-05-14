import matplotlib.pyplot as plt
import argparse as ap
import numpy as np

def load_and_average (filenames):
    summsfA=[]
    nummsfA=0
    for a in filenames:
        thisres,thismsf=np.loadtxt(a,unpack=True)
        if len(summsfA)==0:
            summsfA=thismsf[:]
            nummsfA=1
        else:
            summsfA+=thismsf
            nummsfA+=1
    if nummsfA>0:
        summsfA/=nummsfA
    return thisres,summsfA

p=ap.ArgumentParser()
p.add_argument('-rmsfA',type=str,default=[],action='append',help='set A rmsfs')
p.add_argument('-rmsfB',type=str,default=[],action='append',help='set B rmsfs')
p.add_argument('-rmsfC',type=str,default=[],action='append',help='set C rmsfs')
p.add_argument('-rmsfD',type=str,default=[],action='append',help='set D rmsfs')
p.add_argument('-labelAB',type=str,default='A/B',help='A-B label')
p.add_argument('-labelCD',type=str,default='C/D',help='C-D label')
p.add_argument('-delta-label',type=str,default='',help='delta label')
p.add_argument('-o',type=str,default='plot.png',help='output plot name')

args=p.parse_args()

resA,msfA=load_and_average(args.rmsfA)
resB,msfB=load_and_average(args.rmsfB)
dABmsf=msfA-msfB
resC,msfC=load_and_average(args.rmsfC)
resD,msfD=load_and_average(args.rmsfD)
dCDmsf=msfC-msfD

fig,ax=plt.subplots(1,1,figsize=(6,5))
ax.set_xlabel('Residue Number')
ax.set_ylabel('ΔRMSF {:s} (Å)'.format(args.delta_label))

ax.plot(resA,dABmsf,label=args.labelAB)
ax.plot(resA,dCDmsf,label=args.labelCD)
ax.legend()
plt.savefig(args.o,bbox_inches='tight')
