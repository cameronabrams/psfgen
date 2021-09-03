from mutation import Mutation
from ssbond import SSBond
from graft import Graft
from crot import Crot
from attach import Attach
from link import Link
from deletion import Deletion
from cleavage import Cleavage

ModTypes={'mutations':Mutation,'grafts':Graft,'deletions':Deletion,'crotations':Crot,'attachments':Attach,'links':Link,'ssbonds':SSBond,'cleavages':Cleavage}
comment_chars=['#','%']
class ModsFile:
    def __init__(self,filename=''):
        self.filename=filename
        if len(filename)>0:
            self.allmods={}
            self.delmods={}
            self.activeDict=None
            with open(filename,"r") as fp:
                for l in fp:
                    ll=l.strip()
                    if len(ll)==0:   # blank line; close stanza
                        self.close_stanza()
                    elif ll[0]=='[':  # key: open stanza
                        self.open_stanza(ll)
                    elif ll[0] in comment_chars:
                        pass
                    else:
                        self.process_line(ll)
   
    def open_stanza(self,label=''):
        if len(label)>0:
            if label[1]=='!':
                token=label[2:-1]
                print('Opening negative stanza for {}'.format(token))
                if token in ModTypes:
                    self.current_stanza=token
                    self.delmods[token]=[]
                    self.activeDict=self.delmods
                    self.current_type=ModTypes[token]
            elif label[0]=='[':
                token=label[1:-1]
                print('Opening stanza for {}'.format(token))
                if token in ModTypes:
                    self.current_stanza=token
                    self.allmods[label]=[]
                    self.activeDict=self.allmods
                    self.current_type=ModTypes[label[1:-1]]
#                print('Opening stanza {} with type {}'.format(label,self.current_type))
            else:
                self.current_type=str
                self.activeDict=None
#                print('Opening unknown stanza {}'.format(label))
    def close_stanza(self):
#        print('Closing current stanza {}'.format(self.current_stanza))
        self.current_type='None'
        self.current_stanza=None
    def process_line(self,l):
#        print('Processing line',l,'in',self.current_stanza)
        self.activeDict[self.current_stanza].append(self.current_type(l))
    def show_type(self,typ=None):
        for k in self.allmods.keys():
            stanza=k
#            print('looking for',stanza,'for',typ,'in',ModTypes)
            if stanza in ModTypes:
#                print('Found',ModTypes[stanza])
                if ModTypes[stanza] == typ:
#                    print('returning',k,ModTypes[stanza])
                    return self.allmods[k]
        return []
    def show_excl(self,typ=None):
        for k in self.delmods.keys():
            stanza=k
#            print('looking for',stanza,'for',typ,'in',ModTypes)
            if stanza in ModTypes:
#                print('Found',ModTypes[stanza])
                if ModTypes[stanza] == typ:
#                    print('returning',k,ModTypes[stanza])
                    return self.delmods[k]
        return []
    def report(self):
        print('Mods in {:s}:'.format(self.filename))
        for k,v in self.allmods.items():
            print(k,':',', '.join([str(_) for _ in v]))
        if len(self.delmods)>0:
            print('Delmods in {:s}:'.format(self.filename))
            for k,v in self.delmods.items():
                print(k,':',', '.join([str(_) for _ in v]))
    def write(self):
        retstr=''
        for k,v in self.allmods.items():
            retstr+='['+k+']\n'
            for l in v:
                retstr+=str(l)+'\n'
            retstr+='\n'
        for k,v in self.delmods.items():
            retstr+='!['+k+']\n'
            for l in v:
                retstr+=str(l)+'\n'
            retstr+='\n'
        return retstr
    def __str__(self):
        return self.filename

if __name__=="__main__":
    import sys
    import argparse as ap
    parser=ap.ArgumentParser()
    parser.add_argument('-f',metavar='<name>',help='testing modsfile name')
    parser.add_argument('-n',nargs='+',default=[])
    args=parser.parse_args()
    for k,v in vars(args).items():
        if type(v) is list:
            print('-{:s} '.format(k)+' '.join(v),end=' ')
        else:
            print('-{:s} {}'.format(k,v),end=' ')
    print()
    TryMods=ModsFile(args.f)
    TryMods.report()
    print('show_type',TryMods.show_type(Mutation))
    output=TryMods.write()
    print(output)