#!/usr/bin/env python
# coding: utf-8




import numpy as np
import math
import argparse



parser = argparse.ArgumentParser(description='Window parser 4 ChromoPainter: /n This script combine samples.out and label file to provide window based copyng vectors')
parser.add_argument('--label',help='labelfile used for CP',required=True)
parser.add_argument('--samples',help='samples.out file',required=True)
parser.add_argument('--genmap',help='genetic map in bim format',required=True)
#parser.add_argument('--winlength',help='window length in bases',required=True)
parser.add_argument('--out',help='output name',required=True)

args = vars(parser.parse_args())

label=args['label']
samplefile=args['samples']
mapfile=args['genmap']
outname = args['out']

### EDIT YOUR WINDOW LENGHT HERE
winlength= 500000






def parselabel (label):
    dictlab={}
    with open (label,"r") as l:
        for i,line in enumerate(l):
            dictlab[i+1]=line.strip().split()[0]
    return (dictlab)

    
    





def getPoscM(bim):
    poscM = []
    with open(bim, "r") as b:
        for line in b:
            line = line.strip().split()
            poscM.append(line[2])
    return poscM





def extractsamples(samplefile):
    myinds=[]
    with open (samplefile) as inp:
        for line in inp:
                if line.startswith("HAP"):
                    myinds.append(line.strip().split()[2])
        return(myinds)





def parselabelPops (label):
    dictlab={}
    i=1
    with open (label,"r") as l:
        for line in l:
            line=line.strip().split()
            if line[2]=="1":
                dictlab[i]=line[1]
                dictlab[i+1]=line[1]
            i+=2
    return (dictlab)





def getdonPop(arr,pops):
    arr=map(int,arr)
    arr2=[]
    for i in arr:
        arr2.append(pops[i])
    arr2=np.asarray(arr2)
    return(arr2)





def parseSample(samplefile,poplab):
    arrTot=[]
    with open (samplefile) as inp:
        for line in inp:
                if line.startswith(("HAP","EM_iter")):
                    continue
                else:
                    line=line.strip().split(" ")[1:]
                    arr = np.asarray(line)
                    arrWPop=list(getdonPop(arr, poplab))
                    arrTot.append(arrWPop)
    return(np.asarray(arrTot))





def splithaps(array):
    a=array[0:int(q.shape[0]/2)]
    b=array[int(q.shape[0]/2):]
    return(a,b)





def getPos(bim,winsize):
    pos={}
    with open (bim,"r") as b:
        for i,line in enumerate(b):
            line=line.strip().split()
            pos[i]=line+[math.floor(int(line[3])/winsize)]
    return(pos)





def createwindows(y): 
    allcvs=dict()
    cv=makeEmptyCV(label)
    w=1
    for n in range(1,len(y)-1):
        if pos[n+1][4]==pos[n][4]:
            cv[y[n]]+=1
        elif pos[n+1][4]!=pos[n][4]:
            w+=1
            #print(cv)
            cv=makeEmptyCV(label)
            cv[y[n]]+=1
        allcvs[w]=[x for x in cv.values()]
    return(allcvs)





def makeEmptyCV(label):
    with open (label) as l:
        ll=[]
        for line in l:
            line=line.strip().split()[1]
            if line not in ll:
                ll.append(line)
        return({x:0 for x in ll})





def createwindowscM(y,dist): 
    allcvs=dict()
    cv=makeEmptyCV(label)
    w=1
    for n in range(1,len(y)-1):
        #print(n)
        if pos[n+1][4]==pos[n][4]:
            cv[y[n]]+=dist[n]
        elif pos[n+1][4]!=pos[n][4]:
            w+=1
            #print(cv)
            cv=makeEmptyCV(label)
            cv[y[n]]+=dist[n]
        allcvs[w]=[x for x in cv.values()]
    return(allcvs)





q=parseSample(samplefile,parselabelPops(label))





q=splithaps(q)





pos=getPos(mapfile,winlength)





poscm = getPoscM(mapfile)
dist = [0.0] + [float(poscm[x + 1]) - float(poscm[x]) for x in range(len(poscm) - 1)]





empty=makeEmptyCV(label)




indnames=extractsamples(samplefile)




with open(outname,"w") as o:
    o.write("\t".join(["IND","winstart","winend"]+[x for x in makeEmptyCV(label).keys()])+"\n")
    for nind,inds in enumerate(q):
        i=1
        indhaps=dict()
        for haps in inds:
            #print(haps)
            singlehap=createwindowscM(haps,dist)
            indhaps[i]=singlehap
            i+=1
        for ww in indhaps[1].keys():
            #print(indhaps[2][ww])
            z=[indhaps[x][ww] for x in range(1,len(indhaps))]
            k=sum(map(np.array, z))
            #print(k)
            k=[indnames[nind]] +[str(ww*winlength)]+[str((ww+1)*winlength)]+ list(map(str,k))
            o.write("\t".join(k)+"\n")

