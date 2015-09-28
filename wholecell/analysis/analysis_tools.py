#!/usr/bin/env python

"""
Analysis script toolbox functions
"""

import os
import wholecell.utils.constants

LOW_RES_DIR = 'low_res_plots'
SVG_DIR = 'svg_plots'
LOW_RES_DPI = 120
DEFAULT_IMAGE_TYPE = '.pdf'

def exportFigure(plt, plotOutDir, plotOutFileName, metadata=None):

	if metadata != None:
		if "gen" in metadata:
			# Format metadata signature for single gen figure
			metadata_signature = "_".join([str(metadata["time"])[:14],
					str(metadata["variant_function"]),
					str(metadata["variant_index"]),
					"Gen", str(metadata["gen"])+'/'+str(metadata["total_gens"]),
					"Githash", str(metadata["git_hash"])[:12],
					"Desc", str(metadata["description"])])
		else:
			# Format metadata signature for multi gen figure
			metadata_signature = "_".join([str(metadata["time"][:14]),
					str(metadata["variant_function"]),
					str(metadata["variant_index"]),
					"Total_gens", str(metadata["total_gens"]),
					"Githash", str(metadata["git_hash"])[:12],
					"Desc", str(metadata["description"])])
			
		# Add metadata signature to the bottom of the plot
		plt.figtext(0,0, metadata_signature)

	# Make folders for holding alternate types of images
	if not os.path.exists(os.path.join(plotOutDir, LOW_RES_DIR)):
		os.mkdir(os.path.join(plotOutDir, LOW_RES_DIR))
	if not os.path.exists(os.path.join(plotOutDir, SVG_DIR)):
		os.mkdir(os.path.join(plotOutDir, SVG_DIR))

	# Save PDF image
	plt.savefig(os.path.join(plotOutDir, plotOutFileName + DEFAULT_IMAGE_TYPE))

	# Save SVG image
	plt.savefig(os.path.join(plotOutDir, SVG_DIR, plotOutFileName + '.svg'))

	# Save PNG image
	plt.savefig(os.path.join(plotOutDir, LOW_RES_DIR, plotOutFileName + '.png'), dpi=LOW_RES_DPI)