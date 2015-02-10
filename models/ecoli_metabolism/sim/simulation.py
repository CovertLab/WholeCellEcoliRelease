
import os

from wholecell.sim.simulation import Simulation

# States
from wholecell.states.bulk_molecules import BulkMolecules

# Processes
from models.ecoli_metabolism.processes.transcript_elongation import TranscriptElongation
from models.ecoli_metabolism.processes.rna_degradation import RnaDegradation
from models.ecoli_metabolism.processes.polypeptide_elongation import PolypeptideElongation
from models.ecoli_metabolism.processes.protein_degradation import ProteinDegradation
from models.ecoli_metabolism.processes.replication import Replication
from models.ecoli_metabolism.processes.metabolism import Metabolism
from models.ecoli_metabolism.processes.maintenance import Maintenance

# Listeners
from models.ecoli.listeners.mass import Mass

from models.ecoli_metabolism.listeners.concentration_change import ConcentrationChange

# Analysis
from models.ecoli.analysis.single import massFractions
from models.ecoli.analysis.single import evaluationTime
from models.ecoli.analysis.single import processMassBalance
from models.ecoli_metabolism.analysis.single import effectiveBiomass

# Initialization
from models.ecoli_metabolism.sim.initial_conditions import calcInitialConditions


class EcoliMetabolismSimulation(Simulation):
	_stateClasses = (
		BulkMolecules,
		)

	_processClasses = (
		TranscriptElongation,
		RnaDegradation,
		PolypeptideElongation,
		ProteinDegradation,
		Replication,
		Metabolism,
		Maintenance,
		)

	_listenerClasses = (
		Mass,
		ConcentrationChange
		# NtpUsage, # restore these in some general sense...
		# AAUsage,
		)

	_analysisSingleFiles = (
		massFractions.__file__,
		evaluationTime.__file__,
		processMassBalance.__file__,
		effectiveBiomass.__file__
		)

	_hookClasses = ()

	_initialConditionsFunction = calcInitialConditions

	_lengthSec = 3600
	_timeStepSec = 1

	_logToShell = True
	_shellColumnHeaders = [
		"Time (s)",
		"Dry mass (fg)",
		"Dry mass fold change",
		"Protein fold change",
		"RNA fold change",
		"Expected fold change"
		]

	_logToDisk = False

	@classmethod
	def printAnalysisSingleFiles(cls, fileName = None):
		L = []
		fileList = sorted(cls._analysisSingleFiles)
		if fileName == None:
			for f in fileList:
				if f.endswith(".pyc") or f == "__init__.py":
					print f[:-1]
					L.append(f[:-1])
				else:
					print f
					L.append(f)
		else:
			h = open(fileName, "w")
			for f in fileList:
				if f.endswith(".pyc") or f == "__init__.py":
					h.write(f[:-1] + "\n")
					L.append(f[:-1])
				else:
					h.write(f + "\n")
					L.append(f)
			h.close()

		return L

	@classmethod
	def printAnalysisCohortFiles(cls, fileName = None):
		fileList = []
		if fileName == None:
			for f in fileList:
				if f.endswith(".pyc") or f == "__init__.py":
					print f[:-1]
				else:
					print f
		else:
			h = open(fileName, "w")
			for f in fileList:
				if f.endswith(".pyc") or f == "__init__.py":
					h.write(f[:-1] + "\n")
				else:
					h.write(f + "\n")
			h.close()