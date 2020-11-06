inputFileDir="/home/lu/AcrossTissue/csvs/6LS.csv"
outputFileDIr="/home/lu/AcrossTissue/csvs/6LS_filter.csv"
f=open(inputFileDir,"r+")
lines=f.readlines()
outputString=""
outputF=open(outputFileDIr,"w")
for line in lines:
    if line.split(",").count("")<5:
        outputF.write(line)
outputF.close()
