#!/usr/bin/python

#read in test data (avoid repitition!)

import imp


with open('data/testFile.txt','r') as myfile:
    testData = myfile.read()

splitData = testData.split('Clstr')[1:]


