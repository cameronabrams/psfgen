import argparse
from mutation import Mutation
from molecule import Molecule
mut=[]

parser=argparse.ArgumentParser()
parser.add_argument('pdb',nargs='+',metavar='<?.pdb>',type=Molecule,help='name(s) of pdb file to parse; first is treated as the base molecule')
parser.add_argument('-mut',metavar='X_Y###Z',default=[],action='append',type=Mutation,help='specify mutation.  Format: X is chainID, Y is one-letter residue code to mutate FROM, ### is sequence number, and Z is one-letter residue code to mutate TO.  Multiple -mut\'s can be specified.')
parser.add_argument('-kc',action='store_true',default=False)

args=parser.parse_args()
print(args)

Mol=args.pdb
print(Mol[0])
if len(args.mut)>0:
    for m in args.mut:
        print(m)
if args.kc:
    print('keeping conflicts')
