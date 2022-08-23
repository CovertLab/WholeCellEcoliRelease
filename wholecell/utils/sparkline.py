'''
Utility functions for generating sparkline plots
'''

from __future__ import absolute_import, division, print_function

import numpy as np


def whitePadSparklineAxis(axis, xAxis=True, secondary=False):
	axis.spines["right"].set_visible(secondary)
	axis.spines["top"].set_visible(False)
	axis.spines["left"].set_position(("outward", 10))
	axis.spines["right"].set_position(("outward", 10))
	axis.spines["bottom"].set_position(("outward", 10))
	axis.set_yticks(axis.get_ylim())
	axis.set_xticks(axis.get_xlim())
	axis.tick_params(which="both", direction="out", right=secondary, top=False)

	if not xAxis:
		axis.spines["bottom"].set_visible(False)
		axis.tick_params(bottom = False)
		axis.tick_params(axis = "x", labelbottom=False)

def simpleSparklineAxis(axis):
	axis.spines['top'].set_visible(False)
	axis.spines['bottom'].set_visible(False)
	axis.xaxis.set_ticks_position('none')
	axis.tick_params(which = 'both', direction = 'out')


def sparklineAxis(axis, x, y, tickPos, lineType, color):
	axis.plot(x, y, linestyle=lineType, drawstyle='steps', color=color, linewidth=2)
	axis.spines['top'].set_visible(False)
	axis.spines['bottom'].set_visible(False)
	axis.yaxis.set_ticks_position(tickPos)
	axis.xaxis.set_ticks_position('none')
	axis.tick_params(which = 'both', direction = 'out')

	for tl in axis.get_yticklabels():
		tl.set_color(color)

def setAxisMaxMinY(axis, data):
	if np.isnan(data).all():
		return

	ymax = np.nanmax(data)
	ymin = 0
	if ymin == ymax:
		axis.set_yticks([ymin])
	else:
		axis.set_yticks([ymin, ymax])

def setAxisMaxMinX(axis, data):
	if np.isnan(data).all():
		return

	xmax = np.nanmax(data)
	xmin = 0
	if xmin == xmax:
		axis.set_xticks([xmin])
	else:
		axis.set_xticks([xmin, xmax])
