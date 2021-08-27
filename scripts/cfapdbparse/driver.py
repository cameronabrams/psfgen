import sys
import operator
from datetime import date 
from molecule import Molecule
import argparse
''' 
    Testing parsing of PDB and CIF files
    Cameron F Abrams
    cfa22@drexel.edu
'''

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    print('driver {} / python {}'.format(date.today(),sys.version.replace('\n',' ').split(' ')[0]))
    Molecules=[]

    parser.add_argument('-pdb',action='append',metavar='<?.pdb>',default=[],help='name of PDB file')
    parser.add_argument('-cif',action='append',metavar='<?.cif>',default=[],help='name of CIF file')
    parser.add_argument('-n',nargs='+',default=[])
    args=parser.parse_args()

    PDBfiles=args.pdb
    CIFfiles=args.cif    
    for p in PDBfiles:
        m=Molecule(pdb=p)
        m.summarize()
    for c in CIFfiles:
        m=Molecule(cif=c)
        m.summarize()
