#!/usr/bin/python

#this one is using python3 because it doesn't need anything from dadi, so
#we're free to use the superior python

import sys

if sys.argv[-1] != '--help':
    infile = sys.argv[1]
    with open(infile,'r') as myfile:
        initialData = myfile.read()

    #ok, now we can start to build our new file
    #so dadi wants us to have 

    goodData = initialData.split('Clstr')[1:]

###### dictionary containing population information

populations = {'CMK018': 1,'CMK019': 1,'CMK020': 1,'CMK021': 1,
'KFS165': 1,'KFS166': 1,'KFS167': 1,'KFS168': 1,'KFS169': 1,'KFS170': 1,
'NIR007': 1,'NIR008': 1,'NIR009': 1,'NIR011': 1,'NIR012': 1,'NIR013': 1,
'NIR019': 1,'NIR021': 1,'NIR022': 1,'CMK022': 2,'CMK023': 2,'CMK024': 2,
'CMK025': 2,'KFS171': 2,'KFS172': 2,'KFS174': 2,'KFS175': 2,'KFS176': 2,
'NIR024': 2,'NIR025': 2,'NIR026': 2,'NIR027': 2,'NIR029': 2,'NIR030': 2,
'NIR031': 2,'NIR032': 2,'B29859': 3,'KFS025': 3,'KFS028': 3,'KFS130': 3,
'TPS018': 3,'JSB003': 3,'KFS054': 3,'KFS084': 3,'KFS088': 3,'KFS100': 3,
'MDS125': 3,'KFS195': 3,'CMK015': 3,'KFS142': 3,'KFS186': 3,'KFS164': 3}

###### sorenson data class

class SorensonCluster:
    """
    this class allows for the processing of a single cluster,
    so the way to do it is to iterate through the clusters in
    your file and use the 'outputDadi' method on each cluster
    object
    """
    def __init__(self,data):
        self.lines = data.split('\n')
        self.lines.pop()
        self.snpData = self.lines[1]
        self.stats = self.lines[0]

    def outputDadi(self):
        dadiDict = {1:{},2:{},3:{}}
        alleleDict = {}
        for line in self.lines[2:]:
            splitIt = line.split('\t')
            pop = populations[splitIt[0]]
            alleleNum = splitIt[3]
            if alleleNum in alleleDict:
                pass
            else:
                alleleDict[alleleNum] = splitIt[2]
            if alleleNum in dadiDict[pop]:
                dadiDict[pop][alleleNum] += 1
            else:
                dadiDict[pop][alleleNum] = 1
        prelimAlleles = list(set(
                list(dadiDict[1].keys()) +
                list(dadiDict[2].keys()) +
                list(dadiDict[3].keys())))
        #alleleList = [int(x) for x in prelimAlleles]
        outAllele = max(dadiDict[3], key=dadiDict[3].get)
        inDict = dict(dadiDict[1])
        for key in dadiDict[2].keys():
            if key in inDict:
                inDict[key] += dadiDict[2][key]
            else:
                inDict[key] = dadiDict[2][key]
        inAllele = max(inDict, key=inDict.get)
        return dadiDict, prelimAlleles, alleleDict, outAllele, inAllele
        


            

    


#get some help!

if sys.argv[-1] == '--help':
    print("""Usage:\n
\tmakeSnpFile.py inputFile.out --option\n
Options:
--help\t\tget this menu\n--snpdata\tprint out information about all snps
--stats\t\tget stats for each cluster
--dadi\t\tprint a formatted dadi snp file to stdout (use > file.txt to save to a file)""")


#just for fun, a way to print out data for all snps
if sys.argv[-1] == '--snpdata':
    for i in goodData:
        print(SorensonCluster(i).snpData)

#and for the 'stats' line
if sys.argv[-1] == '--stats':
    for i in goodData:
        print(SorensonCluster(i).stats)


#if sys.argv[2] == 'dadi':





