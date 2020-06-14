#plot each column in a different color in a scatter plot, assume the first line is the header


import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
import scipy.stats as ss
import validator
import matplotlib as mpl
mpl.rcsetup.cycler = cycler(color=[])
colorCycleList=['red', 'green','#e6e600','purple','#3377ff','pink']
colorCycleList=['green','#e6e600','purple','#3377ff','pink']
headers=[]
inputFile="/home/lu/AcrossTissue/ComparisonLSG/selectionParameterBetweenLS.csv"


#format: first line should be the header, each column is a list
def parseInputFile(fname):
    csvMatrix=[]
    f=open(fname)
    lines=f.readlines()
    global headers
    print "Line 0 (Usually the header):",
    print lines[0]
    headers=lines[0].rstrip().split(",")
    for line in lines[1:]:
        splitList=line.rstrip().split(",")
        csvMatrix.append(splitList)
    return csvMatrix



def convertMatirxFromStringToFloat(matrix):
    newMatrix=[]
    for ls in matrix:
        floatList = list(map(float, ls)) 
        newMatrix.append(floatList)
    return newMatrix

csvMatrix=parseInputFile(inputFile)


matrix=np.transpose(np.asarray(csvMatrix))
        
AANames=matrix[0]    
codonNames=matrix[1]
matrix=matrix[2:]
headers=headers[2:]
matrix=convertMatirxFromStringToFloat(matrix)

def median(lst):
    n = len(lst)
    s = sorted(lst)
    return (sum(s[n//2-1:n//2+1])/2.0, s[n//2])[n % 2] if n else None


def getMedianListFromMatrix(matrix,transpose=False):
    medianList=[]
    if transpose:
        matrix=np.transpose(matrix)
    for ls in matrix:
        median=np.median(ls)
        medianList.append(median)
    return medianList

def getMeanListFromMatrix(matrix,transpose=False):
    meanList=[]
    if transpose:
        matrix=np.transpose(matrix)
    for ls in matrix:
        mean=np.mean(ls)
        meanList.append(mean)
    return meanList


#return lst1/lst2, must ensure the length is the same
def divideListByList(lst1,lst2):
    resultList=[]
    for i in range(len(lst1)):
        if lst2[i]==0:
            resultList.append(0)
            continue
        result=lst1[i]/lst2[i]
        resultList.append(result)
    return resultList

        






#works, but the proportion is not clear in visualization
def plotMatrixStackedBar(matrix):
    global colorCycleList
    fig=plt.figure()
    plotGroupLS=[]
    meanList=getMeanListFromMatrix(matrix,transpose=True)
    medianList=getMedianListFromMatrix(matrix,transpose=True)
    
    
    corelation=validator.testCorelation(medianList,meanList,"pearson" )
    if corelation<0.99:
        print "Warning! Mean and Median has a corelation less than 0.99, choice of mean of median might make a difference"
#    plotMean=plt.bar(range(1,len(meanList)+1),meanList)
#    plotGroupLS.append(plotMean)
    
    for i in range(len(matrix)):
        ls=matrix[i]
        proportionLS=divideListByList(ls,meanList)
        #default color map
#        plotGroup=plt.bar(range(1,len(proportionLS)+1),proportionLS,linewidth=0.1,edgecolor='b')
#        #specifiend color map
        pickedColor=colorCycleList[i%len(matrix)]
        plotGroup=plt.bar(range(1,len(proportionLS)+1),proportionLS,color=pickedColor,alpha=0.5,edgecolor=pickedColor,linewidth=0.1)
#        
        
        plotGroupLS.append(plotGroup)
    
    global codonNames
    labelList=[]
    for s in codonNames:
        newS=""
        for c in s:
            newS+=c
            newS+="\n"
        labelList.append(newS)
    plt.xticks(range(1,len(medianList)+1),labelList,rotation=-5, fontsize=5)
    legendGroup=[]
    for i in range(len(plotGroupLS)):
        legendGroup.append(plotGroupLS[i][0])
    plt.legend(legendGroup,headers,fancybox=True,fontsize=5)
    plt.ylabel("Selection Estimate Proportion")
    plt.axhline(y=1, color='grey', linestyle='--',linewidth=0.5)

    plt.show()
    fig.savefig("test.pdf",bbox_inches='tight')

#plotMatrixStackedBar(matrix)

def minusAListByItsAvg(lst):
    avg=np.mean(lst)
    newLst=[]
    for x in lst:
        newLst.append(x-avg)
    return newLst        
    

def plotHistogramByAA(matrix):
    global AANames
    global codonNames
    n=len(headers)
    fig, ax = plt.subplots(n, n,figsize=(40, 24))
    for i in range(len(headers)):
        for j in range(len(headers)):
            
            subAxis=ax[i,j]
            if i<=j:
                subAxis.set_visible(False)
                continue
                
#            print i
#            print headers[i]
#            print j
#            print headers[j]
        #defining sub plots
            
            subAxis.set_title("%s VS %s"%(headers[i],headers[j]))
            LSlist=matrix[i]
            cl=colorCycleList[i%len(colorCycleList)]
            aaPointer=AANames[0]
            labelList=[]
            xCounter=1
            xList=[]
            yList=[]
            
            cnt=0
            for k in range(len(AANames)):
                aa=AANames[k]
                codonName=codonNames[k]
                if aa==aaPointer:      
                    cnt+=1
                else:
#                    print xList
#                    print yList
                    yList=minusAListByItsAvg(yList)
                    subAxis.bar(xList,yList,color=cl,alpha=0.7)
                    xList=[]
                    yList=[]
                    xCounter+=1
                    subAxis.bar(xCounter,0,color=cl)
                    xCounter+=1
                    subAxis.bar(xCounter,0,color=cl)
                    labelList.append("")
                    labelList.append("")

            
                newS=""
                for c in codonName:
                    newS+=c
                    newS+="\n"
                    
                labelList.append(newS)
                yList.append(LSlist[k])
                xCounter+=1
                xList.append(xCounter)
                aaPointer=aa

            
            LSlist=matrix[j]
            cl=colorCycleList[j%len(colorCycleList)]
            aaPointer=AANames[0]
            xCounter=1
            xList=[]
            yList=[]
            for k in range(len(AANames)):
                aa=AANames[k]
                codonName=codonNames[k]
                if aa==aaPointer:      
                    cnt+=1
                else:
#                    print xList
#                    print yList
                    yList=minusAListByItsAvg(yList)
                    subAxis.bar(xList,yList,color=cl,alpha=0.7)
                    xList=[]
                    yList=[]
                    xCounter+=1
                    subAxis.bar(xCounter,0,color=cl)
                    xCounter+=1
                    subAxis.bar(xCounter,0,color=cl)


                yList.append(LSlist[k])
                xCounter+=1
                xList.append(xCounter)
                aaPointer=aa
            subAxis.bar(-1,0,label=headers[i],color=colorCycleList[i%len(colorCycleList)])
            subAxis.bar(-1,0,label=headers[j],color=colorCycleList[j%len(colorCycleList)])
            subAxis.set_xticks(range(len(labelList)))
            subAxis.legend()
            subAxis.set_xticklabels(labelList,rotation=-5, fontsize=5)
            subAxis.set_ylabel("Delta Eta Relative to Mean Across AA Group")

                


                
               


    fig.savefig("test.pdf")       
            
#    plt.bar(xList,yList)
#    plt.xticks(xList,labelList,rotation=-5, fontsize=5)



plotHistogramByAA(matrix)



    
    














