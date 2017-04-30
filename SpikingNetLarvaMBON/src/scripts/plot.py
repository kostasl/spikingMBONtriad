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

# Import PdfPages
from matplotlib.backends.backend_pdf import PdfPages

##Data files with columns:
#t	R	Vm	SpikeRate
list_of_files 	  = [('../../dat/MBONLog.csv','MBON'),('../../dat/DANLog.csv','DAN'),('../../dat/KCLog.csv','KC')]
list_of_filesBase = [('../../dat/MBONLog_Naive.csv','MBON'),('../../dat/DANLog_Naive.csv','DAN'),('../../dat/KCLog_Naive.csv','KC')]

#Report Files being Opened
for filename, label in list_of_files:
	print(filename)

# Use numpy to load the data contained in the file
dataMbon = np.loadtxt(list_of_files[0][0],delimiter='\t',skiprows=1)
dataDAN = np.loadtxt(list_of_files[1][0],delimiter='\t',skiprows=1)
dataKC = np.loadtxt(list_of_files[2][0],delimiter='\t',skiprows=1)

##Select data with R = 0
rates = np.unique(dataMbon[:,1])

# Create plots with pre-defined labels.
fig, ax = plt.subplots()

for R in rates:
        rdat = dataMbon[np.where(dataMbon[:,1] == R)]
	ax.plot(rdat[1::50,0],rdat[1::50,3],'b', label='MBON')
	ax.axis([0, 50, 0, 100])

        rdat = dataDAN[np.where(dataDAN[:,1] == R)]
	ax.plot(rdat[1::50,0],rdat[1::50,3],'r', label='DAN')
	ax.axis([0, 50, 0, 100])

        rdat = dataDAN[np.where(dataDAN[:,1] == R)]
	ax.plot(rdat[1::50,0],rdat[1::50,4],'g', label='R')
	ax.axis([0, 50, 0, 100])

        rdat = dataKC[np.where(dataKC[:,1] == R)]
	ax.plot(rdat[1::50,0],rdat[1::50,3],'k', label='KC')
	ax.axis([0, 50, 0, 100])

plt.xlabel('Time (sec)')
plt.ylabel('Hz')

#legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')

#Get artists and labels for legend and chose which ones to display
handles, labels = ax.get_legend_handles_labels()
display = (0,1,2,3)


#Create legend from custom artist/label lists
ax.legend([handle for i,handle in enumerate(handles) if i in display],
          [label for i,label in enumerate(labels) if i in display])

##Next is for adding Custom Labels and lines
#Create custom artists
#mbonArtist = plt.Line2D((0,1),(0,0), color='b', marker='', linestyle='')
#danArtist = plt.Line2D((0,1),(0,0), color='g')
#KCArtist = plt.Line2D((0,1),(0,0), color='k')
#ax.legend([handle for i,handle in enumerate(handles) if i in display]+[mbonArtist,danArtist,KCArtist],
#          [label for i,label in enumerate(labels) if i in display]+['MBON', 'DAN','KC'])



plt.show()

# Initialize the pdf file
pp = PdfPages('multiFig.pdf')

# Save the figure to the file
pp.savefig(fig)

# Close the file
pp.close()
