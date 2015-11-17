#!/usr/bin/env python

"""
Reusable plotting functions and tools

@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/06/2015
"""

import os
import wholecell.utils.constants

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


def plotSplom(arrayOfdataArrays, nameArray="", plotCorrCoef=True):
	"""
	Plot a scatterplot matrix (Splom) of data contained in arrayOfdataArrays,
	with labels in the same order held within nameArray.
	"""

	if len(arrayOfdataArrays) != len(nameArray):
		raise IndexError("Your array of data arrays and the array of names must be the same length.")

	num_entries = len(arrayOfdataArrays)
	plottingIndex = 1
	for rowNum in xrange(1,num_entries+1):
		for colNum in xrange(1,num_entries+1):
			if colNum < plottingIndex:
				continue
			plt.subplot(num_entries,num_entries,num_entries*(rowNum-1)+(colNum))
			plt.scatter(arrayOfdataArrays[colNum-1], arrayOfdataArrays[rowNum-1])
			
			if nameArray != "":
				plt.xlabel(nameArray[colNum-1])
				plt.ylabel(nameArray[rowNum-1])
			
			if plotCorrCoef:
				corr_coef, pValue = stats.pearsonr(arrayOfdataArrays[colNum-1], arrayOfdataArrays[rowNum-1])
				xLocation = np.amax([.9*np.amax(arrayOfdataArrays[colNum-1]),(1-((np.amax(arrayOfdataArrays[colNum-1]) - np.amin(arrayOfdataArrays[colNum-1]))*.1))*np.amax(arrayOfdataArrays[colNum-1])])
				yLocation = np.amax([.9*np.amax(arrayOfdataArrays[rowNum-1]),(1-((np.amax(arrayOfdataArrays[rowNum-1]) - np.amin(arrayOfdataArrays[rowNum-1]))*.1))*np.amax(arrayOfdataArrays[rowNum-1])])
				plt.text(xLocation,yLocation,"r = %.4f" % (corr_coef), verticalalignment='top', horizontalalignment='left')
		plottingIndex += 1