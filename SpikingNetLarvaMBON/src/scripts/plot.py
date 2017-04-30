#!/usr/bin/python
# coding: utf-8 
#plot Spiking net results
# ploting and datahandling script for SpikingNetLarvaMBON simulation experiments
# this script is run at the end of simulation to plot the Neural Triad responses comparing gain change
# between Naive and trained configurations across different reward inputs.
# The script below selects/sections the Firing rate data according to the R input firing rate
# Author: Kostantinos Lagogiannis 2017

import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

##Data files with columns:
#t	R	Vm	SpikeRate
filename = '../../dat/MBONLog.csv'
filenameBase = '../../dat/MBONLog_Naive.csv'

# Use numpy to load the data contained in the file
dataMbon = np.loadtxt(filename,delimiter='\t',skiprows=1)

##Select data with R = 0
rates = np.unique(dataMbon[:,1])


for R in rates:
        rdat = dataMbon[np.where(dataMbon[:,1] == R)]
	plt.plot(rdat[1::50,0],rdat[1::50,3],'b')
	plt.axis([0, 50, 0, 100])
	plt.xlabel('Time (sec)')
	plt.ylabel('Hz')


plt.show()
