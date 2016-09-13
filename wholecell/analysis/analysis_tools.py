#!/usr/bin/env python

"""
Analysis script toolbox functions
"""

import os
import wholecell.utils.constants

LOW_RES_DIR = 'low_res_plots'
SVG_DIR = 'svg_plots'
HTML_DIR = 'html_plots'
LOW_RES_DPI = 120
DEFAULT_IMAGE_TYPE = '.pdf'

def exportFigure(plt, plotOutDir, plotOutFileName, metadata=None):

	if metadata != None and "analysis_type" in metadata:
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
		elif metadata["analysis_type"] == 'variant':
			# Format metadata signature for variant figure
			metadata_signature = "_".join([str(metadata["time"][:13]),
					str(metadata["total_variants"]), "variants",					
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


import mpld3
from mpld3 import plugins, utils

class MetadataSignature(plugins.PluginBase):
	"""Plugin for getting info on click"""

	JAVASCRIPT = """
	mpld3.register_plugin("metadata_signature", MetadataSignature);
	MetadataSignature.prototype = Object.create(mpld3.Plugin.prototype);
	MetadataSignature.prototype.constructor = MetadataSignature;
	function MetadataSignature(fig, props){
	    mpld3.Plugin.call(this, fig, props);
	};

	MetadataSignature.prototype.draw = function(){
		console.log(this);
	    this.fig.canvas.append("text")
	        .text(this.props.id)
	        .style("font-size", 12)
	        .style("opacity", 0.8)
	        .style("text-anchor", "left")
	        .attr("x", this.fig.width *.11)
	        .attr("y", (this.fig.height * .99))
	}
	"""
	def __init__(self, metadata_string):
		self.dict_ = {"type": "metadata_signature",
						"id": metadata_string}

def exportHtmlFigure(fig, plt, plotOutDir, plotOutFileName, metadata=None):

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
		# plt.figtext(0,0, metadata_signature, figure=fig, verticalalignment='bottom')
		plugins.connect(fig, MetadataSignature(metadata_signature))

	if not os.path.exists(os.path.join(plotOutDir, HTML_DIR)):
		os.mkdir(os.path.join(plotOutDir, HTML_DIR))

	mpld3.save_html(fig, str(os.path.join(plotOutDir, HTML_DIR, plotOutFileName) + '.html'))