#!/usr/bin/python

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

#######
class cluster:
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
        #first count number of alleles and record their positions
        snpLine = self.snpData.split('Indels')[0].split('\t')[1:]
        snpLine.pop()
        numSnps = 0
        clusterNum = self.stats.split('\t')[0].split(' ')[1]
        numSnps = len(snpLine)

        #we loop through snpLine, using the original strategy

        snpStrings = []

        for i in range(len(snpLine)):
            #prelim processing
            curSnp = snpLine[i]
            if curSnp == 'none':
                break
            if ',' in curSnp:
                useSnp = curSnp.split(',')[0]
            else:
                useSnp = curSnp

            #count snp states
            dadiDict = {1:{}, 2:{}, 3:{}}
            alleleList = []

            for line in self.lines[2:]:
                splitIt = line.split('\t')
                if splitIt[1] == '.':
                    continue

                pop = populations[splitIt[0]]

                #allele counts
                snpStates = splitIt[2]
                if snpStates[i] in dadiDict[pop]:
                    dadiDict[pop][snpStates[i]] += 1
                else:
                    dadiDict[pop][snpStates[i]] = 1

                #allele states
                if snpStates[i] not in alleleList:
                    alleleList.append(snpStates[i])

            #figure out the stuff we need to print
            
            outAllele = max(dadiDict[3], key=dadiDict[3].get)
            inDict = dict(dadiDict[1])
            for key in dadiDict[2].keys():
                if key in inDict:
                    inDict[key] += dadiDict[2][key]
                else:
                    inDict[key] = dadiDict[2][key]
            inAllele = max(inDict, key=inDict.get)

            #ensure we'll report 0 allele freqs. instead of nothing
            for pop in [1,2,3]:
                for allele in alleleList:
                    if allele not in dadiDict[pop]:
                        dadiDict[pop][allele] = 0

            snpStrings.append(str(inAllele + '\t' + outAllele + '\t' + \
                    alleleList[0] + '\t' + \
                    str(dadiDict[1][alleleList[0]]) + '\t' + \
                    str(dadiDict[2][alleleList[0]]) + '\t' + \
                    str(dadiDict[3][alleleList[0]]) + '\t' + \
                    alleleList[1] + '\t' + \
                    str(dadiDict[1][alleleList[1]]) + '\t' + \
                    str(dadiDict[2][alleleList[1]]) + '\t' + \
                    str(dadiDict[3][alleleList[1]]) + '\t' + \
                    clusterNum + '\t' + useSnp))
        
        finalSnpString = '\n'.join(snpStrings)
        if finalSnpString == '':
            pass
        else:
            return finalSnpString


