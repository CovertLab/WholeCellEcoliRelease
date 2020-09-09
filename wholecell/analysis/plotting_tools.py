"""
Reusable plotting functions and tools
"""

from __future__ import absolute_import, division, print_function

import matplotlib.pyplot as plt
from scipy import stats
from six.moves import range


COLORS_LARGE = ["#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
		"#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
		"#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
		"#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
		"#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
		"#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
		"#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
		"#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",

		"#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
		"#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
		"#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
		"#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
		"#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
		"#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
		"#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
		"#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58"]

COLORS_SMALL = ["#FF0000", "#00FF00", "#0000FF", "#FF00FF", "#00FFFF", "#000000", "#007FFF",
		"#236B8E", "#70DB93", "#B5A642", "#5F9F9F", "#B87333", "#2F4F2F", "#9932CD", "#871F78", "#855E42",
		"#545454", "#8E2323", "#238E23", "#CD7F32", "#527F76",
		"#9F9F5F", "#8E236B", "#2F2F4F", "#CFB53B", "#FF7F00", "#DB70DB",
		"#5959AB", "#8C1717", "#238E68", "#6B4226", "#8E6B23", "#00FF7F",
		"#38B0DE", "#DB9370", "#5C4033", "#4F2F4F", "#CC3299", "#99CC32"]

CMAP_COLORS_255 = [
	[247,247,247],
	[209,229,240],
	[146,197,222],
	[67,147,195],
	[33,102,172],
	[5,48,97],
	]

COLORS_256 = [ # From colorbrewer2.org, qualitative 8-class set 1
	[228,26,28],
	[55,126,184],
	[77,175,74],
	[152,78,163],
	[255,127,0],
	[255,255,51],
	[166,86,40],
	[247,129,191]
	]

def plotSplom(arrayOfdataArrays, nameArray="", stdArrays=None, labels=None, fig=None, plotCorrCoef=True, formatString='o'):
	"""
	Plot a scatterplot matrix (Splom) of data contained in arrayOfdataArrays,
	with labels in the same order held within nameArray.
	"""

	if len(arrayOfdataArrays) != len(nameArray):
		raise IndexError("Your array of data arrays and the array of names must be the same length.")

	if stdArrays is None:
		stdArrays = [None]*len(arrayOfdataArrays)

	if len(stdArrays) != len(arrayOfdataArrays):
		raise IndexError("If you provide an array of standard deviations, there must be one entry per input data array. Entries can be None.")

	if fig is None:
		fig = plt.figure()

	num_entries = len(arrayOfdataArrays)
	plottingIndex = 1
	for rowNum in range(1,num_entries+1):
		for colNum in range(1,num_entries+1):
			if colNum < plottingIndex:
				continue
			plt.subplot(num_entries,num_entries,num_entries*(rowNum-1)+(colNum))

			plt.errorbar(arrayOfdataArrays[colNum-1], arrayOfdataArrays[rowNum-1], xerr=stdArrays[colNum-1], yerr=stdArrays[rowNum-1], fmt=formatString)

			if nameArray != "":
				plt.xlabel(nameArray[colNum-1])
				plt.ylabel(nameArray[rowNum-1])

			if plotCorrCoef:
				corr_coef, pValue = stats.pearsonr(arrayOfdataArrays[colNum-1], arrayOfdataArrays[rowNum-1])
				plt.title("R = %.4f" % (corr_coef,))
		plottingIndex += 1

	return fig
