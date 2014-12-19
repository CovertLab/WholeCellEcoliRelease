
import os

from wholecell.sim.simulation import Simulation

# States
from wholecell.states.bulk_molecules import BulkMolecules
from wholecell.states.unique_molecules import UniqueMolecules
from wholecell.states.bulk_chromosome import BulkChromosome

# Processes
from models.ecoli.processes.complexation import Complexation
from models.ecoli.processes.metabolism import Metabolism
from models.ecoli.processes.rna_degradation import RnaDegradation
from models.ecoli.processes.replication import Replication
from models.ecoli.processes.polypeptide_initiation import PolypeptideInitiation
from models.ecoli.processes.polypeptide_elongation import PolypeptideElongation
from models.ecoli.processes.transcript_initiation import TranscriptInitiation
from models.ecoli.processes.transcript_elongation import TranscriptElongation
from models.ecoli.processes.protein_degradation import ProteinDegradation
from models.ecoli.processes.atp_usage import AtpUsage

# Listeners
from models.ecoli.listeners.mass import Mass
from models.ecoli.listeners.replication_fork import ReplicationForkPosition
from models.ecoli.listeners.ntp_usage import NtpUsage
from models.ecoli.listeners.aa_usage import AAUsage
from models.ecoli.listeners.ribosome_data import RibosomeData
from models.ecoli.listeners.gene_copy_number import GeneCopyNumber
from models.ecoli.listeners.unique_molecule_counts import UniqueMoleculeCounts
from models.ecoli.listeners.fba_results import FBAResults

# Analysis
import models.ecoli.analysis.single
import models.ecoli.analysis.cohort

from models.ecoli.sim.initial_conditions import calcInitialConditions

class EcoliSimulation(Simulation):
	_stateClasses = (
		BulkMolecules,
		UniqueMolecules,
		BulkChromosome
		)

	_processClasses = (
		Metabolism,
		RnaDegradation,
		TranscriptInitiation,
		TranscriptElongation,
		PolypeptideInitiation,
		PolypeptideElongation,
		Replication,
		ProteinDegradation,
		Complexation,
		AtpUsage
		)

	_listenerClasses = (
		Mass,
		ReplicationForkPosition,
		NtpUsage,
		AAUsage,
		RibosomeData,
		GeneCopyNumber,
		UniqueMoleculeCounts,
		FBAResults
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
		directory = os.path.dirname(models.ecoli.analysis.single.__file__)
		fileList = sorted(os.listdir(directory))
		if fileName == None:
			for f in fileList:
				if f.endswith(".pyc") or f == "__init__.py":
					continue
				print os.path.join(directory, f)
		else:
			h = open(fileName, "w")
			for f in fileList:
				if f.endswith(".pyc") or f == "__init__.py":
					continue
				h.write(os.path.join(directory, f) + "\n")
			h.close()

	@classmethod
	def printAnalysisCohortFiles(cls, fileName = None):
		directory = os.path.dirname(models.ecoli.analysis.cohort.__file__)
		fileList = sorted(os.listdir(directory))
		if fileName == None:
			for f in fileList:
				if f.endswith(".pyc") or f == "__init__.py":
					continue
				print os.path.join(directory, f)
		else:
			h = open(fileName, "w")
			for f in fileList:
				if f.endswith(".pyc") or f == "__init__.py":
					continue
				h.write(os.path.join(directory, f) + "\n")
			h.close()