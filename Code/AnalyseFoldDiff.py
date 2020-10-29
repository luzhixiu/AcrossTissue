

import readCSVfiles as RCF
import numpy as np
import matplotlib.pyplot as plt



matrix,rowHeaderList,columnHeaderList= RCF.readCSV("/home/lu/AcrossTissue/Code/LifeStageGenes.csv",readAsString=True)


def floatify(lst):
    return [float(x) for x in lst]
    
def sortX_by_Y(X,Y):
    return [x for _,x in sorted(zip(Y,X))]

#be careful these lists are string typed, for numeric analysis, remember to convert them to float
LS_List=matrix[0]
LS_ExpList=floatify(matrix[1])
secondMaxList=floatify(matrix[2])
restMeanList=floatify(matrix[3])
foldDiffList=floatify(matrix[4])
geneIdList=columnHeaderList
from cycler import cycler

import matplotlib as mpl

colorCycle=(plt.rcParams['axes.prop_cycle'].by_key()['color'])


colorDict=dict()
def syncColors(typeList):
    global colorDict
    colors=colorCycle[:len(typeList)]
    for x in range(len(typeList)):
        colorDict[typeList[x]]=colors[x]
    return colorDict
typeList=list(set(LS_List))

colorDict=syncColors(typeList)

#get the percentage of each life stage
def getPecentage(LS_List, show_pie=True,titleText=""):
    global colorDict
    typeList=list(set(LS_List))
    percentList=[]
    for ls in typeList:
        cnt=LS_List.count(ls)
        per=cnt/float(len(LS_List))
        percentList.append(per)


    
    colours=[]

    percentList=sorted(percentList)
    typeList=sortX_by_Y(typeList,percentList)
    for x in typeList:
        colours.append(colorDict[x])
        
    import matplotlib.pyplot as plt

    labels = typeList
    plt.figure(dpi=150)
    patches, texts = plt.pie(percentList, startangle=90,colors=colours)
    plt.legend(patches, labels,loc="upper right")
    plt.axis('equal')
    plt.title(titleText)
    plt.tight_layout()
    plt.show()

    return  percentList,typeList

#star visualizing the data a bit
getPecentage(LS_List,titleText="Life Stage Gene Pool (Size 6922)")

#Cut off of foldDiff at 10
foldDiffCutOff=10
cuttedlsLst=[]
for i in range(len(LS_List)):
    if foldDiffList[i]>=foldDiffCutOff:
        cuttedlsLst.append(LS_List[i])
print (len(cuttedlsLst))
getPecentage(cuttedlsLst,titleText="Life Stage Gene Pool (FoldDiffCutOff=10, Size=2977)")

#top 1000
getPecentage(LS_List[:1000],titleText="Life Stage Gene Pool (Top 1000)")

#top 500
getPecentage(LS_List[:500],titleText="Life Stage Gene Pool (Top 500)")

#top 300
getPecentage(LS_List[:300],titleText="Life Stage Gene Pool (Top 300)")

#top 100
getPecentage(LS_List[:100],titleText="Life Stage Gene Pool (Top 100)")

#top 50
getPecentage(LS_List[:30],titleText="Life Stage Gene Pool (Top 50)")








matrix=np.transpose(matrix)


class geneObj():
    def __init__(self, GeneId,LS,LS_EXP,secondMax,restMean,foldDiff):
        self.GeneId = GeneId
        self.LS = LS
        self.LS_EXP=LS_EXP
        self.secondMax=secondMax
        self.restMean=restMean
        self.foldDiff=foldDiff
        
geneObjDict=dict()        
for i in range(len(geneIdList)):
    geneId=geneIdList[i]
    ls=LS_List[i]
    lsExp=float(LS_ExpList[i])
    secondMax=float(secondMaxList[i])
    restMean=float(restMeanList[i])
    foldDiff=float(foldDiffList[i])
    geneId=geneIdList[i]
    geneObjDict[geneIdList[i]]=geneObj(geneId,ls,lsExp,secondMax,restMean,foldDiff)



print()


    
    
        
        