from __future__ import absolute_import


import os
import re
import cPickle
import time

import numpy as np

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.utils import parallelization, units

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from scipy.stats import pearsonr

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS
from models.ecoli.analysis import variantAnalysisPlot

SHUFFLE_VARIANT_TAG = "ShuffleParams"
PLACE_HOLDER = -1

FONT_SIZE=9
trim = 0.05


def getPCCProteome((variant, ap, monomerIds, schmidt_counts)):
	try:
		simDir = ap.get_cells(variant = [variant])[0]
		simOutDir = os.path.join(simDir, "simOut")
		sim_data = cPickle.load(open(ap.get_variant_kb(variant), "rb"))

		ids_translation = sim_data.process.translation.monomerData["id"].tolist()
		schmidt_idx = [ids_translation.index(x) for x in monomerIds]

		monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
		avgCounts = monomerCounts.readColumn("monomerCounts").mean(axis=0)
		sim_schmidt_counts = avgCounts[schmidt_idx]

		pcc, pval = pearsonr(np.log10(sim_schmidt_counts + 1), np.log10(schmidt_counts + 1))

		return pcc, pval

	except Exception as e:
		print e
		return np.nan, np.nan


def getPCCFluxome((variant, ap, toyaReactions, toyaFluxesDict, toyaStdevDict)):

	try:
		simDir = ap.get_cells(variant = [variant])[0]

		sim_data = cPickle.load(open(ap.get_variant_kb(variant), "rb"))
		cellDensity = sim_data.constants.cellDensity

		simOutDir = os.path.join(simDir, "simOut")

		massListener = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = massListener.readColumn("cellMass")
		dryMass = massListener.readColumn("dryMass")
		massListener.close()

		coefficient = dryMass / cellMass * sim_data.constants.cellDensity.asNumber(MASS_UNITS / VOLUME_UNITS)

		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
		reactionFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (fbaResults.readColumn("reactionFluxes").T / coefficient).T
		fbaResults.close()

		modelFluxes = {}
		for toyaReaction in toyaReactions:
			fluxTimeCourse = []

			for rxn in reactionIDs:
				if re.findall(toyaReaction, rxn):
					reverse = 1
					if re.findall("(reverse)", rxn):
						reverse = -1

					if len(fluxTimeCourse):
						fluxTimeCourse += reverse * reactionFluxes[:, np.where(reactionIDs == rxn)]
					else:
						fluxTimeCourse = reverse * reactionFluxes[:, np.where(reactionIDs == rxn)]

			if len(fluxTimeCourse):
				if toyaReaction not in modelFluxes:
					modelFluxes[toyaReaction] = []
				modelFluxes[toyaReaction].append(np.mean(fluxTimeCourse).asNumber(units.mmol / units.g / units.h))

		toyaVsReactionAve = []
		for rxn, toyaFlux in toyaFluxesDict.iteritems():
			if rxn in ["ISOCITDEH-RXN", "SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31."]:
				continue
			if rxn in modelFluxes:
				toyaVsReactionAve.append((np.mean(modelFluxes[rxn]), toyaFlux.asNumber(units.mmol / units.g / units.h), np.std(modelFluxes[rxn]), toyaStdevDict[rxn].asNumber(units.mmol / units.g / units.h)))

		toyaVsReactionAve = np.array(toyaVsReactionAve)
		pcc, pval = pearsonr(toyaVsReactionAve[:,0], toyaVsReactionAve[:,1])

		return pcc, pval
	except Exception as e:
		print e
		return np.nan, np.nan

def getDivisionTime((variant, ap)):
	try:
		simDir = ap.get_cells(variant = [variant])[0]

		simOutDir = os.path.join(simDir, "simOut")

		time_column = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")

		return (time_column.max() - initialTime) / 60.
	except Exception as e:
		print e
		return np.nan

def getInitialMass((variant, ap)):
	try:
		simDir = ap.get_cells(variant = [variant])[0]

		simOutDir = os.path.join(simDir, "simOut")

		mass = TableReader(os.path.join(simOutDir, "Mass"))
		cellDry = mass.readColumn("dryMass")
		return cellDry[0]
	except Exception as e:
		print e
		return np.nan

def getFinalMass((variant, ap)):
	try:
		simDir = ap.get_cells(variant = [variant])[0]

		simOutDir = os.path.join(simDir, "simOut")

		mass = TableReader(os.path.join(simOutDir, "Mass"))
		cellDry = mass.readColumn("dryMass")
		return cellDry[-1]
	except Exception as e:
		print e
		return np.nan


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if metadata is not None and SHUFFLE_VARIANT_TAG not in metadata["variant"]:
			print "This plot only runs for variants where parameters are shuffled."
			return

		if not os.path.isdir(inputDir):
			raise Exception, "variantDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		validation_data = cPickle.load(open(validationDataFile, "rb"))

		schmidtCounts = validation_data.protein.schmidt2015Data["glucoseCounts"]

		toyaReactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
		toyaFluxes = validation_data.reactionFlux.toya2010fluxes["reactionFlux"]
		toyaStdev = validation_data.reactionFlux.toya2010fluxes["reactionFluxStdev"]
		toyaFluxesDict = dict(zip(toyaReactions, toyaFluxes))
		toyaStdevDict = dict(zip(toyaReactions, toyaStdev))

		ap = AnalysisPaths(inputDir, variant_plot = True)

		pool = parallelization.pool(processes=self.cpus)

		# Get simulation time data
		start = time.time()
		args = zip(range(ap.n_variant), [ap] * ap.n_variant)
		divisionResult = pool.map(getDivisionTime, args)
		stop = time.time()
		cPickle.dump(divisionResult, open(os.path.join(plotOutDir, plotOutFileName + "_division.cPickle"), "w"))
		print "%d seconds:\tTo get simulation time data (to compute division time) -- completed" % (stop - start)

		# Get initial mass
		start = time.time()
		args = zip(range(ap.n_variant), [ap] * ap.n_variant)
		initialMassResult = pool.map(getInitialMass, args)
		stop = time.time()
		cPickle.dump(initialMassResult, open(os.path.join(plotOutDir, plotOutFileName + "_initialMass.cPickle"), "w"))
		print "%d seconds:\tTo get initial mass data -- completed" % (stop - start)

		# Get final mass
		start = time.time()
		args = zip(range(ap.n_variant), [ap] * ap.n_variant)
		finalMassResult = pool.map(getFinalMass, args)
		stop = time.time()
		cPickle.dump(finalMassResult, open(os.path.join(plotOutDir, plotOutFileName + "_finalMass.cPickle"), "w"))
		print "%d seconds:\tTo get final mass data -- completed" % (stop - start)

		# Get fluxome correlation
		start = time.time()
		args = zip(range(ap.n_variant), [ap] * ap.n_variant, [toyaReactions] * ap.n_variant, [toyaFluxesDict] * ap.n_variant, [toyaStdevDict] * ap.n_variant)
		fluxomeResult = pool.map(getPCCFluxome, args)
		stop = time.time()
		cPickle.dump(fluxomeResult, open(os.path.join(plotOutDir, plotOutFileName + "_fluxome.cPickle"), "w"))
		print "%d seconds:\tTo get fluxome correlation -- completed" % (stop - start)

		# Get proteome correlation
		start = time.time()
		args = zip(range(ap.n_variant), [ap] * ap.n_variant, [validation_data.protein.schmidt2015Data["monomerId"].tolist()] * ap.n_variant, [schmidtCounts] * ap.n_variant)
		proteomeResult = pool.map(getPCCProteome, args)
		stop = time.time()
		cPickle.dump(proteomeResult, open(os.path.join(plotOutDir, plotOutFileName + "_proteome.cPickle"), "w"))
		print "%d seconds:\tTo get proteome correlation -- completed" % (stop - start)

		pool.close()
		pool.join()

		divisionResult = np.array(divisionResult)
		initialMassResult = np.array(initialMassResult)
		finalMassResult = np.array(finalMassResult)
		fluxomeResult = np.array(fluxomeResult)
		proteomeResult = np.array(proteomeResult)

		dataTable = np.hstack((divisionResult.reshape(-1, 1), proteomeResult, fluxomeResult, initialMassResult.reshape(-1, 1), finalMassResult.reshape(-1, 1)))
		colNames = ["division time", "proteome pearson r", "proteome pearson r p value", "fluxome pearson r", "fluxome pearson r p value", "initial dry mass (fg)", "final dry mass (fg)"]
		h = open(os.path.join(plotOutDir, "distribution_division_fluxome_proteome_data_matrix.tsv"), "w")

		h.write("\t".join(colNames) + "\n")
		np.savetxt(h, dataTable, delimiter = "\t")
		h.close()


if __name__ == "__main__":
	Plot().cli()
