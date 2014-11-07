#!/usr/bin/python

#this one is using python3 because it doesn't need anything from dadi, so
#we're free to use the superior python

import sys
from cluster import cluster

###### read in data
if sys.argv[-1] != '--help':
    infile = sys.argv[1]
    with open(infile,'r') as myfile:
        initialData = myfile.read()

    # split the data by cluster
    goodData = initialData.split('Clstr')[1:]

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
        print(cluster(i).snpData)

#and for the 'stats' line
if sys.argv[-1] == '--stats':
    for i in goodData:
        print(cluster(i).stats)


if sys.argv[2] == 'dadi':
    for i in goodData:
        print(cluster(i).outputDadi())





