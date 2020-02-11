"""
Analysis script toolbox functions
"""

import os

import numpy as np

from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath

LOW_RES_DIR = 'low_res_plots'
SVG_DIR = 'svg_plots'
HTML_DIR = 'html_plots'
LOW_RES_DPI = 120

def exportFigure(plt, plotOutDir, plotOutFileName, metadata=None, transparent=False,
        dpi=LOW_RES_DPI, extension=None):

	if metadata is not None and "analysis_type" in metadata:
		analysis_type = metadata["analysis_type"]

		if analysis_type == 'single':
			# Format metadata signature for single gen figure
			metadata_signature = "_".join([
				str(metadata["time"])[:13],
				str(metadata["variant_function"]),
				str(metadata["variant_index"]),
				"Seed",
				str(metadata["seed"]),
				"Gen",
				str(metadata["gen"]) + '/' + str(int(metadata["total_gens"]) - 1),
				"Githash",
				str(metadata["git_hash"])[:10],
				"Desc",
				str(metadata["description"])
				])

		elif analysis_type == 'multigen':
			# Format metadata signature for multi gen figure
			metadata_signature = "_".join([
				str(metadata["time"][:13]),
				str(metadata["variant_function"]),
				str(metadata["variant_index"]),
				"Seed",
				str(metadata["seed"]),
				str(metadata["total_gens"]),
				"gens",
				"Githash",
				str(metadata["git_hash"])[:10],
				"Desc",
				str(metadata["description"])
				])

		elif analysis_type == 'cohort':
			# Format metadata signature for cohort figure
			metadata_signature = "_".join([
				str(metadata["time"][:13]),
				str(metadata["variant_function"]),
				str(metadata["variant_index"]),
				str(metadata["total_gens"]),
				"gens",
				"Githash",
				str(metadata["git_hash"])[:10],
				"Desc",
				str(metadata["description"])
				])

		elif analysis_type == 'variant':
			# Format metadata signature for variant figure
			metadata_signature = "_".join([
				str(metadata["time"][:13]),
				str(metadata["total_variants"]),
				"variants",
				str(metadata["total_gens"]),
				"gens",
				"Githash",
				str(metadata["git_hash"])[:10],
				"Desc",
				str(metadata["description"])
				])

		else:
			raise ValueError('Unknown analysis_type {}'.format(analysis_type))

		# Add metadata signature to the bottom of the plot
		# Don't accidentally trigger $TeX formatting$.
		metadata_signature = metadata_signature.replace('$', '')
		plt.figtext(0,0, metadata_signature, size=8)

	# Make folders for holding alternate types of images
	filepath.makedirs(plotOutDir, LOW_RES_DIR)
	filepath.makedirs(plotOutDir, SVG_DIR)

	# Save images
	if extension:
		# Only save one type in main analysis directory if extension is given
		plt.savefig(os.path.join(plotOutDir, plotOutFileName + extension), dpi=dpi, transparent=transparent)
	else:
		# Save all image types
		plt.savefig(os.path.join(plotOutDir, plotOutFileName + '.pdf'), transparent=transparent)
		plt.savefig(os.path.join(plotOutDir, SVG_DIR, plotOutFileName + '.svg'), transparent=transparent)
		plt.savefig(os.path.join(plotOutDir, LOW_RES_DIR, plotOutFileName + '.png'), dpi=dpi, transparent=transparent)

def read_bulk_molecule_counts(sim_out_dir, mol_names):
	'''
	Reads a subset of molecule counts from BulkMolecules using the indexing method
	of readColumn. Should only be called once per simulation being analyzed with
	all molecules of interest.

	Args:
		sim_out_dir (str): path to the directory with simulation output data
		mol_names (list-like or tuple of list-like): lists of strings containing
			names of molecules to read the counts for. A single array will be
			converted to a tuple for processing.

	Returns:
		generator of ndarray: int counts with all time points on the first dimension
		and each molecule of interest on the second dimension. The number of
		generated arrays will be separated based on the input dimensions of mol_names
		(ie if mol_names is a tuple of two arrays, two arrays will be generated).

	Example use cases:
		names1 = ['ATP[c]', 'AMP[c]']
		names2 = ['WATER[c]']

		# Read one set of molecules
		(counts1,) = read_bulk_molecule_counts(sim_out_dir, names1)

		# Read two or more sets of molecules
		(counts1, counts2) = read_bulk_molecule_counts(sim_out_dir, (names1, names2))

	TODO: generalize to any TableReader, not just BulkMolecules, if readColumn method
	is used for those tables.
	'''

	# Convert an array to tuple to ensure correct dimensions
	if not isinstance(mol_names, tuple):
		mol_names = (mol_names,)

	# Check for string instead of array since it will cause mol_indices lookup to fail
	for names in mol_names:
		if isinstance(names, basestring):
			raise Exception('mol_names must be a tuple of arrays not strings like {}'.format(names))

	bulk_reader = TableReader(os.path.join(sim_out_dir, 'BulkMolecules'))

	bulk_molecule_names = bulk_reader.readAttribute("objectNames")
	mol_indices = {mol: i for i, mol in enumerate(bulk_molecule_names)}

	lengths = [len(names) for names in mol_names]
	indices = np.hstack([[mol_indices[mol] for mol in names] for names in mol_names])
	bulk_counts = bulk_reader.readColumn2D('counts', indices)

	start_slice = 0
	for length in lengths:
		counts = bulk_counts[:, start_slice:start_slice + length].squeeze()
		start_slice += length
		yield counts
