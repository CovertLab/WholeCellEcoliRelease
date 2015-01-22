#!/usr/bin/env python

import os
import ast
import argparse
import subprocess
import datetime
import cPickle

from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS
import reconstruction.ecoli.knowledge_base
import reconstruction.ecoli.fitkb1
from models.ecoli_metabolism.sim.simulation import EcoliMetabolismSimulation
import wholecell.utils.constants

from models.ecoli.analysis.single import massFractions
from models.ecoli.analysis.single import evaluationTime
from models.ecoli.analysis.single import processMassBalance
from models.ecoli_metabolism.analysis.single import effectiveBiomass

def submissionTime():
	now = datetime.datetime.now()
	return "%s%02d%02d.%02d%02d%02d.%d" % (now.year, now.month, now.day, now.hour, now.minute, now.second, now.microsecond)

def writeStringToFile(fileName, string):
	h = open(fileName, "w")
	h.write(string)
	h.close()

def saveMetadata(baseDir, description):
	metadataDir = os.path.join(baseDir, "metadata")
	os.makedirs(metadataDir)

	gitHash = subprocess.Popen(["git", "rev-parse", "HEAD"], stdout = subprocess.PIPE).communicate()[0]
	gitBranch = subprocess.Popen(["git", "symbolic-ref", "--short", "HEAD"], stdout = subprocess.PIPE).communicate()[0]
	gitDiff = subprocess.Popen(["git", "diff"], stdout = subprocess.PIPE).communicate()[0]

	print "Saving metadata"
	writeStringToFile(os.path.join(metadataDir, "git_hash"), gitHash)
	writeStringToFile(os.path.join(metadataDir, "git_branch"), gitBranch)
	writeStringToFile(os.path.join(metadataDir, "git_diff"), gitDiff)
	writeStringToFile(os.path.join(metadataDir, "description"), description + "\n")

def saveKb(baseDir):
	print "Instantiating unfit knowledgebase"
	kb = reconstruction.ecoli.knowledge_base.KnowledgeBaseEcoli()

	print "Fitting knowledgebase (Level 1)"
	reconstruction.ecoli.fitkb1.fitKb_1(kb)

	kbDirectory = os.path.join(baseDir, "kb")
	os.makedirs(kbDirectory)

	fileName = wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
	kbLocation = os.path.join(kbDirectory, fileName)

	cPickle.dump(kb, open(kbLocation, "wb"), protocol = cPickle.HIGHEST_PROTOCOL)

	return kbLocation


def main(**kwargs):

	baseDir = os.path.join("out", submissionTime())

	saveMetadata(baseDir, kwargs["description"])
	del kwargs["description"]

	kbLocation = saveKb(baseDir)

	kwargs["outputDir"] = os.path.join(baseDir, "%06d" % kwargs["seed"], "simOut")
	kwargs["logToDisk"] = True
	kwargs["overwriteExistingFiles"] = False
	kwargs["kbLocation"] = kbLocation

	os.makedirs(kwargs["outputDir"])
	
	sim = EcoliMetabolismSimulation(**kwargs)

	sim.run()

	plotDir = os.path.join(baseDir, "%06d" % kwargs["seed"], "plotOut")

	massFractions.main(kwargs["outputDir"], plotDir, "massFractions", kbLocation)
	evaluationTime.main(kwargs["outputDir"], plotDir, "evaluationTime", kbLocation)
	processMassBalance.main(kwargs["outputDir"], plotDir, "processMassBalance", kbLocation)
	effectiveBiomass.main(kwargs["outputDir"], plotDir, "effectiveBiomass", kbLocation)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	for key in sorted(DEFAULT_SIMULATION_KWARGS):
		argType = type(DEFAULT_SIMULATION_KWARGS[key])

		if argType == type(None):
			argType = str
		if argType == bool:
			argType = ast.literal_eval

		parser.add_argument(
			"--%s" % key,
			default = DEFAULT_SIMULATION_KWARGS[key],
			type = argType
			)

	parser.add_argument("--description", default = "", type = str)

	args = parser.parse_args().__dict__

	main(**args)
