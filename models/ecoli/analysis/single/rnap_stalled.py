"""
Plot count of stalled RNA polymerases from stalled transcription elongation
"""

import os
from matplotlib import pyplot as plt
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		rnap_data_reader = TableReader(os.path.join(simOutDir, 'RnapData'))
		stalled_rnap_counts = rnap_data_reader.readColumn('didStall')

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time

		plt.plot(time / 60., stalled_rnap_counts)
		plt.xlabel("Time (min)")
		plt.ylabel("Protein Counts")
		plt.title("Stalled RNA Polymerase")
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
