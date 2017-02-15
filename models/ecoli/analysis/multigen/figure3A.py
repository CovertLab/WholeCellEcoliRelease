#!/usr/bin/env python
"""
Plot fluxes for metabolic map figure during a shift

@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/13/17
"""

import argparse
import os
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis

START = 8000
SHIFT = 11000
END = 14000
BURNIN = 15

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get all cells
	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
	allDir = ap.get_cells()

	sim_data = cPickle.load(open(simDataFile, "rb"))
	rxnStoich = sim_data.process.metabolism.reactionStoich

	reactants = [
		"GLC-6-P[c]",
		"FRUCTOSE-6P[c]",
		"FRUCTOSE-16-DIPHOSPHATE[c]",
		"DIHYDROXY-ACETONE-PHOSPHATE[c]",
		"GAP[c]",
		"DPG[c]",
		"G3P[c]",
		"2-PG[c]",
		"PHOSPHO-ENOL-PYRUVATE[c]",
		"PYRUVATE[c]",
		"ACETYL-COA[c]",
		"CIT[c]",
		"CIS-ACONITATE[c]",
		"THREO-DS-ISO-CITRATE[c]",
		"2-KETOGLUTARATE[c]",
		"SUC-COA[c]",
		"SUC[c]",
		"FUM[c]",
		"MAL[c]",
		"GLC-6-P[c]",
		"D-6-P-GLUCONO-DELTA-LACTONE[c]",
		"CPD-2961[c]",
		"RIBULOSE-5P[c]",
		"RIBULOSE-5P[c]",
		"XYLULOSE-5-PHOSPHATE[c]",
		"D-SEDOHEPTULOSE-7-P[c]",
		"FRUCTOSE-6P[c]",
		"3-P-SERINE[c]",
		"SER[c]",
		"ACETYLSERINE[c]",
		"HOMO-CYS[c]",
		"GLT[c]",
		"GLT[c]",
		"SER[c]",
		"HISTIDINAL[c]",
		"PYRUVATE[c]",
		"GLT[c]",
		"GLT[c]",
		"GLT[c]",
		"GLT[c]",
		"L-ASPARTATE[c]",
		"O-PHOSPHO-L-HOMOSERINE[c]",
		"MESO-DIAMINOPIMELATE[c]",
		"L-ARGININO-SUCCINATE[c]",
		"2-KETOGLUTARATE[c]",
		"GLT[c]",
		"L-DELTA1-PYRROLINE_5-CARBOXYLATE[c]",
		"DGMP[c]",
		"DGDP[c]",
		"DGMP[c]",
		"DAMP[c]",
		"DADP[c]",
		"DAMP[c]",
		"TMP[c]",
		"TDP[c]",
		"TMP[c]",
		"DUMP[c]",
		"DUDP[c]",
		"DUMP[c]",
		"DCMP[c]",
		"DCDP[c]",
		"DCMP[c]",
		"GMP[c]",
		"GDP[c]",
		"GMP[c]",
		"AMP[c]",
		"ADP[c]",
		"AMP[c]",
		"UMP[c]",
		"UDP[c]",
		"UMP[c]",
		"CMP[c]",
		"CDP[c]",
		"CMP[c]",
		"GDP[c]",
		"GTP[c]",
		"ADP[c]",
		"ATP[c]",
		"TMP[c]",
		"UDP[c]",
		"UTP[c]",
		"UTP[c]",
		"CDP[c]",
		"CTP[c]",
		"GUANINE[c]",
		"ADENINE[c]",
		"URACIL[c]",
		]

	products = [
		"FRUCTOSE-6P[c]",
		"FRUCTOSE-16-DIPHOSPHATE[c]",
		"DIHYDROXY-ACETONE-PHOSPHATE[c]",
		"GAP[c]",
		"DPG[c]",
		"G3P[c]",
		"2-PG[c]",
		"PHOSPHO-ENOL-PYRUVATE[c]",
		"PYRUVATE[c]",
		"ACETYL-COA[c]",
		"CIT[c]",
		"CIS-ACONITATE[c]",
		"THREO-DS-ISO-CITRATE[c]",
		"2-KETOGLUTARATE[c]",
		"SUC-COA[c]",
		"SUC[c]",
		"FUM[c]",
		"MAL[c]",
		"OXALACETIC_ACID[c]",
		"D-6-P-GLUCONO-DELTA-LACTONE[c]",
		"CPD-2961[c]",
		"RIBULOSE-5P[c]",
		"XYLULOSE-5-PHOSPHATE[c]",
		"RIBOSE-5P[c]",
		"D-SEDOHEPTULOSE-7-P[c]",
		"ERYTHROSE-4P[c]",
		"ERYTHROSE-4P[c]",
		"SER[c]",
		"GLY[c]",
		"CYS[c]",
		"MET[c]",
		"TYR[c]",
		"PHE[c]",
		"TRP[c]",
		"HIS[c]",
		"L-ALPHA-ALANINE[c]",
		"VAL[c]",
		"LEU[c]",
		"ILE[c]",
		"L-ASPARTATE[c]",
		"ASN[c]",
		"THR[c]",
		"LYS[c]",
		"ARG[c]",
		"GLT[c]",
		"GLN[c]",
		"PRO[c]",
		"DGDP[c]",
		"DGTP[c]",
		"DGTP[c]",
		"DADP[c]",
		"DATP[c]",
		"DATP[c]",
		"TDP[c]",
		"TTP[c]",
		"TTP[c]",
		"DUDP[c]",
		"DUTP[c]",
		"DUTP[c]",
		"DCDP[c]",
		"DCTP[c]",
		"DCTP[c]",
		"GDP[c]",
		"GTP[c]",
		"GTP[c]",
		"ADP[c]",
		"ATP[c]",
		"ATP[c]",
		"UDP[c]",
		"UTP[c]",
		"UTP[c]",
		"CDP[c]",
		"CTP[c]",
		"CTP[c]",
		"DGDP[c]",
		"DGTP[c]",
		"DADP[c]",
		"DATP[c]",
		"DUMP[c]",
		"DUDP[c]",
		"DUTP[c]",
		"CTP[c]",
		"DCDP[c]",
		"DCTP[c]",
		"GMP[c]",
		"AMP[c]",
		"UMP[c]",
		]

	plt.figure(figsize = (17, 22))

	for simDir in allDir:
		simOutDir = os.path.join(simDir, "simOut")

		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		if initialTime > END or time[-1] < START:
			continue

		fbaReader = TableReader(os.path.join(simOutDir, "FBAResults"))
		flux = fbaReader.readColumn("reactionFluxes")
		reactionIDs = fbaReader.readAttribute("reactionIDs")
		fbaReader.close()

		for idx, (reactant, product) in enumerate(zip(reactants, products)):
			ax = plt.subplot(10, 9, idx + 1)
			totalFlux = np.zeros_like(flux[:, 0])

			for rxn in rxnStoich:
				if reactant in rxnStoich[rxn] and product in rxnStoich[rxn]:
					if rxnStoich[rxn][reactant] < 0 and rxnStoich[rxn][product] > 0:
						direction = 1
					elif rxnStoich[rxn][reactant] > 0 and rxnStoich[rxn][product] < 0:
						direction = -1
					else:
						continue

					if rxn in reactionIDs:
						totalFlux += flux[:, reactionIDs.index(rxn)] * direction
						if rxn + " (reverse)" in reactionIDs:
							totalFlux -= flux[:, reactionIDs.index(rxn + " (reverse)")] * direction

			timeIdx = np.logical_and(np.logical_or(np.logical_and(time >= START, time < SHIFT), np.logical_and(time > SHIFT + BURNIN, time <= END)), time > initialTime + BURNIN)
			
			ax.axhline(0, color = "#aaaaaa")
			ax.axvline(SHIFT, color = "#aaaaaa")
			ax.plot(time[timeIdx], totalFlux[timeIdx], color = "k")
			ax.set_title("%s to %s" % (reactant, product), fontsize = 4)
			ax.tick_params(axis = "both", labelsize = 4)
			ax.set_axis_off()


	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")

if __name__ == "__main__":
	defaultSimDataFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--simDataFile", help = "KB file name", type = str, default = defaultSimDataFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
