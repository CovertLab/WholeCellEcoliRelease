"""
Comparison of model-derived surface area calculations and molecule
count-derived surface area calculations of the inner and outer leaflets
of the outer membrane.
"""

from __future__ import absolute_import, division, print_function

from six.moves import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader
from wholecell.analysis.plotting_tools import CMAP_COLORS_255


CMAP_COLORS = [[shade/255. for shade in color] for color in CMAP_COLORS_255]
CMAP_OVER = [0, 1, 0.75]

'''
References:
- width of 0.5 um: Neidhart et al., Physiology of the Bacterial Cell, Chapter 1
	- strain: B/R
- width of 0.73 um: Cluzel et al., Nucleic Acids Research (2008)
	- strain: K-12 Frag1
- surface areas per molecule: Harvard Bionumbers (Properties of the major 
components of outer membrane)
	-sources: Neidhart et al., Escherichia coli and Salmonella: Cellular and 
	Molecular Biology, Chapter 8 ; Smit et al., Journal of Bacteriology (1975)
- assumption that half of all phospholipids produced end up in outer membrane:
Osborn et al., Journal of Biological Chemistry (1972)
- surface area validation point of 6 um^2: Neidhart et al., Physiology of the 
Bacterial Cell, Chapter 1
	- strain: B/R  
- average e. coli time point 44% of the way through cell cycle: Neidhart et al., 
Physiology of the Bacterial Cell, Chapter 1

Note:
- cell volume: this value from the model is derived from the density parameter
	- strain: ML308
'''

WIDTH = [0.5, 0.73]		# in um
AVERAGE_SURFACE_AREA = 6	# in um^2

SURFACE_AREA_PER_MOLECULE = {		# in um^2 / molecule
	'LPS': 1.42E-6,
	'porins_and_ompA': 9E-6,
	'phospholipids': 4.71E-07,
	'lipoprotein': 7.14E-07,
}

OUTER_MEM_PROTEIN = {
	# phospholipids
	'CPD-12819[c]': 'phosphatidylethanolamine',
    'CPD-12824[c]' : 'cardiolipin',
    'CPD-8260[c]': 'phosphatidylglycerol',

	# LPS
    'CPD0-939[c]': 'LPS',

	# porins and ompA
    'EG10669-MONOMER[i]': 'ompA',
	'CPLX0-7533[e]': 'outer membrane porin C',
	'CPLX0-7534[e]': 'outer membrane porin F',

	# lipoprotein
	'EG10544-MONOMER[o]': 'murein lipoprotein',
}


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)
		average_timepoint = np.log(
			sim_data.mass.avg_cell_to_initial_cell_conversion_factor) / np.log(2)
		outer_mem_protein_ids = list(OUTER_MEM_PROTEIN.keys())

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		mass = TableReader(os.path.join(simOutDir, "Mass"))

		# Load data
		volume = mass.readColumn("cellVolume")
		initial_time = main_reader.readAttribute('initialTime')
		time = (main_reader.readColumn('time') - initial_time) / 60
		(counts,) = read_bulk_molecule_counts(simOutDir, (outer_mem_protein_ids,))
		counts = counts.astype(float).T


		# Calculate surface area based off of model volume
		radii = np.divide(WIDTH, 2)
		surface_area_model = np.zeros((np.shape(radii)[0], np.shape(volume)[0]))
		for i, radius in enumerate(radii):
			surface_area_sphere = 4 * np.pi * radius**2
			length = (volume - ((4 / 3) * np.pi * radius**3))/(np.pi * radius**2) + 2*radius
			surface_area_cylinder = 2 * np.pi * radius * (length - 2*radius)
			surface_area_model[i, :] = surface_area_cylinder + surface_area_sphere

		# calculate SA based off of molecule counts
		surface_area_LPS = SURFACE_AREA_PER_MOLECULE['LPS'] * counts[outer_mem_protein_ids.index(
			'CPD0-939[c]')]
		surface_area_porins_and_ompA = SURFACE_AREA_PER_MOLECULE['porins_and_ompA'] * np.sum(
			counts[[outer_mem_protein_ids.index('CPLX0-7533[e]'),
					outer_mem_protein_ids.index('CPLX0-7534[e]'),
					outer_mem_protein_ids.index('EG10669-MONOMER[i]')], :], axis = 0)
		surface_area_phospholipids = SURFACE_AREA_PER_MOLECULE['phospholipids'] * 0.5 * np.sum(
			counts[[outer_mem_protein_ids.index('CPD-12819[c]'),
					outer_mem_protein_ids.index('CPD-12824[c]'),
					outer_mem_protein_ids.index('CPD-8260[c]')], :], axis = 0)
		surface_area_lipoprotein = SURFACE_AREA_PER_MOLECULE['lipoprotein'] * counts[
			outer_mem_protein_ids.index('EG10544-MONOMER[o]')]

		surface_area_outer_leaflet = np.add(surface_area_LPS,
			surface_area_porins_and_ompA)
		surface_area_inner_leaflet = np.add(np.add(surface_area_phospholipids,
			surface_area_lipoprotein), surface_area_porins_and_ompA)

		plt.figure(figsize=(8.5, 11))
		plt.plot(time, surface_area_model[0], 'r-')
		plt.plot(time, surface_area_model[1], 'r--')
		plt.plot(time, surface_area_outer_leaflet, 'b-')
		plt.plot(time, surface_area_inner_leaflet, 'b--')
		plt.scatter(average_timepoint * time[-1], AVERAGE_SURFACE_AREA)

		plt.xlabel('time (min)')
		plt.ylabel(r'surface area ($\mu m ^2$)')
		plt.title('Comparison of surface area calculations')
		model_derived_legend = ['Model derived SA: width = ' + str(width)
								for width in WIDTH]
		plt.legend(model_derived_legend + [
			'Molecule count derived SA of outer leaflet',
			'Molecule count derived SA of inner leaflet',
			'Average E. coli SA'])

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
