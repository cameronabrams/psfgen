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
p.add_argument('-rmsfA',type=str,default=[],nargs='+',help='set A rmsfs')
p.add_argument('-rmsfB',type=str,default=[],nargs='+',help='set B rmsfs')
p.add_argument('-rmsfC',type=str,default=[],nargs='+',help='set C rmsfs')
p.add_argument('-rmsfD',type=str,default=[],nargs='+',help='set D rmsfs')
p.add_argument('-resid-wts',type=str,default='',help='file of resid weights')
p.add_argument('-D',type=str,default='A-B',help='delta operation direction')
p.add_argument('-labelAB',type=str,default='A/B',help='A-B label')
p.add_argument('-labelCD',type=str,default='C/D',help='C-D label')
p.add_argument('-delta-label',type=str,default='',help='delta label')
p.add_argument('-o',type=str,default='plot.png',help='output plot name')

args=p.parse_args()
postname,prename=['post','pre']
if '-' in args.delta_label:
    postname,prename = args.delta_label.split('-')

resA,msfA=load_and_average(args.rmsfA)
resB,msfB=load_and_average(args.rmsfB)

if args.D == 'A-B':
    dABmsf=msfA-msfB
else:
    dABmsf=msfB-msfA
resC,msfC=load_and_average(args.rmsfC)
resD,msfD=load_and_average(args.rmsfD)
if args.D == 'A-B':
    dCDmsf=msfC-msfD
else:
    dCDmsf=msfD-msfC


fig,ax=plt.subplots(3,1,figsize=(6,10),sharex=True)
ax[2].set_xlabel('Residue Number')
ax[0].set_ylabel('ΔRMSF {:s} (Å)'.format(args.delta_label))
ax[1].set_ylabel('RMSF (Å)')
ax[2].set_ylabel('RMSF (Å)')

ssAB=dABmsf.mean()
ssCD=dCDmsf.mean()
if len(args.resid_wts)>0:
    wtr,wtv=np.loadtxt(args.resid_wts,unpack=True)
    ax[0].bar(wtr,-wtv/2,bottom=0,label="epitope occupancy",color='black',width=1.0,alpha=0.5)
    ax[0].bar(wtr,wtv/2,bottom=0,label=None,color='black',width=1.0,alpha=0.5)
    ssAB=0.0
    ssCD=0.0
    count=0
    for r,v in zip(wtr,wtv):
        ssAB+=dABmsf[np.flatnonzero(resA==r)[0]]*v
        ssCD+=dCDmsf[np.flatnonzero(resA==r)[0]]*v
        count+=1
    ssAB/=count
    ssCD/=count
# to do: compute the epitope-occupancy-weighted sum-of-fluctuations
ax[0].plot(resA,dABmsf,label=args.labelAB+' ({:.2f} Å$^2$)'.format(ssAB))
ax[0].plot(resA,dCDmsf,label=args.labelCD+' ({:.2f} Å$^2$)'.format(ssCD))
ax[0].legend()
ax[0].set_ylim([-3,3])
ax[1].set_ylim([0,4])
ax[1].plot(resA,msfA,label=(postname if args.D=='A-B' else prename)+'-'+args.labelAB)
ax[1].plot(resA,msfB,label=(prename if args.D=='A-B' else postname)+'-'+args.labelAB)
ax[1].legend()
ax[2].set_ylim([0,4])
ax[2].plot(resA,msfC,label=(postname if args.D=='A-B' else prename)+'-'+args.labelCD)
ax[2].plot(resA,msfD,label=(prename if args.D=='A-B' else postname)+'-'+args.labelCD)
ax[2].legend()
plt.savefig(args.o,bbox_inches='tight')
