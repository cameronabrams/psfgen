def CIFMakeStructs(db):
    structs={}
    for k in db.keys():
        kk=k.split('.')
        s=kk[0]
        if s in structs:
            i+=1
            structs[s][kk[1]]=i
        else:
            i=0
            structs[s]={}
            structs[s][kk[1]]=i
    return structs

def GetCIFStructDictList(db,structs,struk):
    keys=[_ for _ in structs[struk].keys()]
    dlist=[]
    chek=struk+'.'+keys[0]
    isloop=True if db.FindLoop(chek)!=-1 else False
    if isloop:
        x=db.GetLoop(chek)
        for y in x:
            d={}
            for v in keys:
                d[v]=y[structs[struk][v]]
            dlist.append(d)
    else:
        d={}
        for v in keys:
            d[v]=db[struk+'.'+v]
        dlist.append(d)
    return dlist