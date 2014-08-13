
from wholecell.sim.simulation import simulationFactory

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
from models.ecoli.listeners.ribosome_stalling import RibosomeStalling
from models.ecoli.listeners.gene_copy_number import GeneCopyNumber
from models.ecoli.listeners.unique_molecule_counts import UniqueMoleculeCounts
from wholecell.listeners.evaluation_time import EvaluationTime
from models.ecoli.listeners.fba_results import FBAResults

from models.ecoli.sim.initial_conditions import calcInitialConditions

EcoliSimulation = simulationFactory(
	states = [BulkMolecules, UniqueMolecules, BulkChromosome],
	processes = [
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
	], 
	listeners = [
		Mass,
		ReplicationForkPosition,
		NtpUsage,
		AAUsage,
		RibosomeStalling,
		GeneCopyNumber,
		UniqueMoleculeCounts,
		EvaluationTime,
		FBAResults
	],
	hooks = [], # same as default
	initialConditionsFunction = calcInitialConditions,
	lengthSec = 3600, # same as default
	timeStepSec = 1, # same as default
	logToShell = True, # same as default
	shellColumnHeaders = [
		"Time (s)",
		"Dry mass (fg)",
		"Dry mass fold change",
		"Protein fold change",
		"RNA fold change",
		"Expected fold change"
	],
	logToDisk = False, # same as default
	)

sim = EcoliSimulation() # basic sim

sim = EcoliSimulation(seed = 1, lengthSec = 100)

sim = EcoliSimulation(logToDisk = True, outputDir = ..., logToDiskEvery = 10)


"""

Factory (required):

states
processes
listeners
kbLocation
initialConditionsFunction

Factory (optional):

hooks (default ())
lengthSec (default 3600)
timeStepSec (default 1)
logToShell (default True)
shellColumnHeaders (default ?)
logToDisk (default False)
outputDir (default None)
logToDiskEvery (default 1)
overwriteExistingFiles (default False)

Instantiation (optional):

kbLocation
seed
lengthSec
logToShell
logToDisk
outputDir
logToDiskEvery
overwriteExistingFiles

"""
