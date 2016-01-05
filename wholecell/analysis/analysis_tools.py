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
		if metadata["analysis_type"] == 'single':
			# Format metadata signature for single gen figure
			metadata_signature = "_".join([str(metadata["time"])[:13],
					str(metadata["variant_function"]),
					str(metadata["variant_index"]),
					"Seed", str(metadata["seed"]),
					"Gen", str(metadata["gen"])+'/'+str(int(metadata["total_gens"])-1),
					"Githash", str(metadata["git_hash"])[:10],
					"Desc", str(metadata["description"])])
		elif metadata["analysis_type"] == 'multigen':
			# Format metadata signature for multi gen figure
			metadata_signature = "_".join([str(metadata["time"][:13]),
					str(metadata["variant_function"]),
					str(metadata["variant_index"]),
					"Seed", str(metadata["seed"]),
					str(metadata["total_gens"]), "gens",
					"Githash", str(metadata["git_hash"])[:10],
					"Desc", str(metadata["description"])])
		elif metadata["analysis_type"] == 'cohort':
			# Format metadata signature for cohort figure
			metadata_signature = "_".join([str(metadata["time"][:13]),
					str(metadata["variant_function"]),
					str(metadata["variant_index"]),
					str(metadata["total_gens"]), "gens",
					"Githash", str(metadata["git_hash"])[:10],
					"Desc", str(metadata["description"])])

		# Add metadata signature to the bottom of the plot
		plt.figtext(0,0, metadata_signature, size=8)

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