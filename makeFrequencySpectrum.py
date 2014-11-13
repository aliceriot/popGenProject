#!/usr/bin/python2

#note that the #! line is for my (ben's) machine, you may need to change it
#if you're on mac OS or windows

import dadi

#we make a dictionary of the snp data
dd = dadi.Misc.make_data_dict('hunsteinAsnpFile.txt')

#now we can make a spectrum with this data
fs = dadi.Spectrum.from_data_dict(dd, ['1','2','3'],[20,20,20])

#make one we can plot
newfs = dadi.Spectrum.from_data_dict(dd, ['1','2'],[20,20])

import pylab
dadi.Plotting.plot_single_2d_sfs(newfs, vmin=0.1)
pylab.show()
