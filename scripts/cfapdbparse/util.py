def vmd_instructions(fp,script,logname='tmp.log',args='',msg=''):
    fp.write('echo "VMD) script={} log={} msg: {}"\n'.format(script,logname,msg))
    if args!='':
        fp.write(r'$VMD -dispdev text -e '+script+r' -args '+args+r' > '+logname+' 2>&1\n')
    else:
        fp.write(r'$VMD -dispdev text -e '+script+r' > '+logname+' 2>&1\n')
    fp.write('if [ $? -ne 0 ]; then\n')
    fp.write('   echo "VMD failed.  Check the log file {}. Exiting."\n'.format(logname))
    fp.write('   exit 1\n')
    fp.write('fi\n')

def namd_instructions(fp,cfgname,psf,coor,outname,logname,
                      npe=8,numminsteps=0,numsteps=0,seed=0,template='vac.namd',
                      temperature=310,extras=[],msg='',stdparamfiles=[],localparamfiles=[],
                      stdcharmmdir=r'\$env(HOME)/charmm/toppar',
                      localcharmmdir=r'\$env(PSFGEN_BASEDIR)/charmm'):
    fp.write('cat $PSFGEN_BASEDIR/templates/{}'.format(template))
    fp.write('  | sed s/%OUT%/{}/g'.format(outname))
    fp.write('  | sed s/%NUMMIN%/{}/'.format(numminsteps))
    fp.write('  | sed s/%NUMSTEPS%/{}/'.format(numsteps))
    fp.write('  | sed s/%SEED%/{}/g'.format(seed))
    fp.write('  | sed s/%TEMPERATURE%/{}/g'.format(temperature))
    fp.write('  | sed "/#### SYSTEM CONFIGURATION FILES END/i structure {}"'.format(psf))
    fp.write('  | sed "/#### SYSTEM CONFIGURATION FILES END/i coordinates {}"'.format(coor))
    sentinelline='#### PARAMETER FILES END'
    for st in stdparamfiles:
        fp.write(' | sed "/{}/i parameters {}/{}" '.format(sentinelline,stdcharmmdir,st))
    for st in localparamfiles:
        fp.write(' | sed "/{}/i parameters {}/{}" '.format(sentinelline,localcharmmdir,st))
    sentinelline='#### EXTRAS END'
    for ex in extras:
        fp.write('  | sed "/'+sentinelline+'/i '+ex+'" ')
    fp.write(' > {}\n'.format(cfgname))
    namdp='+p{:d}'.format(npe)
    fp.write('echo "NAMD2) config={} log={} outputname={} msg={}"\n'.format(cfgname,logname,outname,msg))
    fp.write(r'$CHARMRUN '+namdp+r' $NAMD2 '+cfgname+r' > '+logname+'\n')
    fp.write('if [ $? -ne 0 ]; then\n')
    fp.write('   echo "NAMD failed.  Check log file {}. Exiting."\n'.format(logname))
    fp.write('   exit 1\n')
    fp.write('fi\n')
    
def MrgCmdLineAndFileContents(cl_list,filename,typ):
    if filename!='':
        with open(filename,'r') as f:
           for l in f:
               if l[0]!='#':
                   cl_list.append(typ(l))
    return cl_list

def DictFromString(string):
    #print('parsing {}'.format(string))
    my_dict = {}
    if len(string)>0:
        items=string.split(',')
        for i in items:
            kv=i.split('=')
            k=kv[0]
            v=kv[1]
            my_dict[k]=v
    return my_dict

def DefOrDict(d,varname,default):
    return default if varname not in d else d[varname]
