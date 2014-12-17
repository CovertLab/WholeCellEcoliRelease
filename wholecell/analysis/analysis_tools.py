#!/usr/bin/env python

'''
Analysis script toolbox functions
'''

import os
import wholecell.utils.constants

LOW_RES_DIR = 'low_res_plots'
PDF_DIR = 'pdf_plots'
LOW_RES_DPI = 30
DEFAULT_IMAGE_TYPE = '.svg'

def exportFigure(plt, plotOutDir, plotOutFileName):

	if not os.path.exists(os.path.join(plotOutDir, LOW_RES_DIR)):
		os.mkdir(os.path.join(plotOutDir, LOW_RES_DIR))
	if not os.path.exists(os.path.join(plotOutDir, PDF_DIR)):
		os.mkdir(os.path.join(plotOutDir, PDF_DIR))

	# Save SVG image
	plt.savefig(os.path.join(plotOutDir, plotOutFileName + DEFAULT_IMAGE_TYPE))

	# Save PDF image
	plt.savefig(os.path.join(plotOutDir, PDF_DIR, plotOutFileName + '.pdf'))

	# Save PNG image
	plt.savefig(os.path.join(plotOutDir, LOW_RES_DIR, plotOutFileName + '.png'), dpi=LOW_RES_DPI)