"""

Defines a class for storing and managing the parameters involved in simulation
instantiation.

"""

from __future__ import division

import collections
import os
import wholecell.utils.constants

import json

# References to sub-simulation abstractions

# States
import wholecell.states.bulk_molecules
import wholecell.states.unique_molecules
# import wholecell.states.chromosome
# import wholecell.states.transcripts
import wholecell.states.bulk_chromosome

STATE_CLASSES = (
	wholecell.states.bulk_molecules.BulkMolecules,
	wholecell.states.unique_molecules.UniqueMolecules,
	# wholecell.states.chromosome.Chromosome,
	# wholecell.states.transcripts.Transcripts,
	wholecell.states.bulk_chromosome.BulkChromosome
	)

STATES = {stateClass.name():stateClass for stateClass in STATE_CLASSES}

# Processes
import wholecell.processes.complexation
import wholecell.processes.metabolism
import wholecell.processes.metabolism_fba
import wholecell.processes.rna_degradation
import wholecell.processes.replication
import wholecell.processes.translation.unique_polypeptide_initiation
import wholecell.processes.translation.unique_polypeptide_elongation
import wholecell.processes.transcription.unique_transcript_initiation
import wholecell.processes.transcription.unique_transcript_elongation
import wholecell.processes.protein_degradation
import wholecell.processes.replication_initiation
import wholecell.processes.biomass_internment
import wholecell.processes.atp_usage

PROCESS_CLASSES = (
	wholecell.processes.complexation.Complexation,
	wholecell.processes.metabolism.Metabolism,
	wholecell.processes.metabolism_fba.MetabolismFba,
	wholecell.processes.rna_degradation.RnaDegradation,
	wholecell.processes.replication.Replication,
	wholecell.processes.translation.unique_polypeptide_initiation.UniquePolypeptideInitiation,
	wholecell.processes.translation.unique_polypeptide_elongation.UniquePolypeptideElongation,
	wholecell.processes.transcription.unique_transcript_initiation.UniqueTranscriptInitiation,
	wholecell.processes.transcription.unique_transcript_elongation.UniqueTranscriptElongation,
	wholecell.processes.protein_degradation.ProteinDegradation,
	wholecell.processes.replication_initiation.ReplicationInitiation,
	wholecell.processes.biomass_internment.BiomassInternment,
	wholecell.processes.atp_usage.AtpUsage,
	)

PROCESSES = {processClass.name():processClass for processClass in PROCESS_CLASSES}

# Listeners
import wholecell.listeners.mass
import wholecell.listeners.metabolic_flux
import wholecell.listeners.replication_fork
import wholecell.listeners.ntp_usage
import wholecell.listeners.aa_usage
import wholecell.listeners.ribosome_stalling
import wholecell.listeners.gene_copy_number
import wholecell.listeners.metabolic_demands
import wholecell.listeners.unique_molecule_counts
import wholecell.listeners.evaluation_time
import wholecell.listeners.effective_biomass_objective

LISTENER_CLASSES = (
	wholecell.listeners.mass.Mass,
	wholecell.listeners.metabolic_flux.MetabolicFlux,
	wholecell.listeners.replication_fork.ReplicationForkPosition,
	wholecell.listeners.ntp_usage.NtpUsage,
	wholecell.listeners.aa_usage.AAUsage,
	wholecell.listeners.ribosome_stalling.RibosomeStalling,
	wholecell.listeners.gene_copy_number.GeneCopyNumber,
	wholecell.listeners.metabolic_demands.MetabolicDemands,
	wholecell.listeners.unique_molecule_counts.UniqueMoleculeCounts,
	wholecell.listeners.evaluation_time.EvaluationTime,
	wholecell.listeners.effective_biomass_objective.EffectiveBiomassObjective
	)

LISTENERS = {listenerClass.name():listenerClass for listenerClass in LISTENER_CLASSES}

# Loggers
import wholecell.loggers.shell
import wholecell.loggers.disk

# TODO: logger logic more consistent with listeners/states/processes

DEFAULT_SHELL_COLUMN_HEADERS = [
	"Time (s)",
	"Dry mass (fg)",
	"Dry mass fold change",
	"Protein fold change",
	"RNA fold change",
	"Expected fold change"
	]

# Hooks
import wholecell.hooks.rnap_count_hook
import wholecell.hooks.ribosome_count_hook

HOOK_CLASSES = (
	wholecell.hooks.rnap_count_hook.RnapCountHook,
	wholecell.hooks.ribosome_count_hook.RibosomeCountHook
	)

HOOKS = {hookClass.name():hookClass for hookClass in HOOK_CLASSES}

# Default parameters

DEFAULT_STATES = (
	'BulkMolecules',
	'UniqueMolecules',
	'BulkChromosome'
	)

DEFAULT_PROCESSES = (
	'Metabolism',
	'RnaDegradation',
	'UniqueTranscriptInitiation',
	'UniqueTranscriptElongation',
	'UniquePolypeptideInitiation',
	'UniquePolypeptideElongation',
	'Replication',
	'ProteinDegradation',
	'Complexation',
	# 'BiomassInternment',
	'AtpUsage'
	)

DEFAULT_LISTENERS = (
	'Mass',
	'ReplicationForkPosition',
	'NtpUsage',
	'AAUsage',
	'RibosomeStalling',
	'GeneCopyNumber',
	'UniqueMoleculeCounts',
	'EvaluationTime',
	# 'EffectiveBiomassObjective'
	)

DEFAULT_HOOKS = ( # NOTE: there should probably never be any default hooks
	)

DEFAULT_LENGTH = 3600 # sec
DEFAULT_TIME_STEP = 1 # sec

DEFAULT_SEED = 0

DEFAULT_KB_LOCATION = os.path.join(
	wholecell.utils.constants.SERIALIZED_KB_DIR,
	wholecell.utils.constants.SERIALIZED_KB_FIT_FILENAME
	)

# TODO: restore KB reconstruction option when available
# TODO: define defaults for reconstruction, logging

# Simulation keywords

SIM_KWARG_DEFAULTS = dict(
	states = DEFAULT_STATES,
	processes = DEFAULT_PROCESSES,
	listeners = DEFAULT_LISTENERS,
	hooks = DEFAULT_HOOKS,
	lengthSec = DEFAULT_LENGTH, timeStepSec = DEFAULT_TIME_STEP,
	seed = DEFAULT_SEED,
	logToShell = True, shellColumnHeaders = DEFAULT_SHELL_COLUMN_HEADERS,
	logToDisk = False, outputDir = None, overwriteExistingFiles = False, logToDiskEvery = None,
	kbLocation = DEFAULT_KB_LOCATION
	)

# Exceptions

class SimDefinitionException(Exception):
	pass


# Classes

class SimDefinition(object):
	"""
	A class for providing convenient access to the static definition of a 
	simulation's parameters, and to perform some of the simpler parsing of 
	parameters into usable objects.
	"""

	def __init__(self, **kwargs):
		# Generic attribute setting
		attributeNames = SIM_KWARG_DEFAULTS.viewkeys()

		for attribute, value in kwargs.iteritems():
			if attribute not in attributeNames:
				raise SimDefinitionException(
					"Unrecognized simulation definition attribute: {}".format(attribute)
					)

			else:
				setattr(self, attribute, value)

		for unassignedAttribute in (attributeNames - kwargs.viewkeys()):
			setattr(self, unassignedAttribute, SIM_KWARG_DEFAULTS[unassignedAttribute])


	def createStates(self):
		return collections.OrderedDict([
			(className, STATES[className]())
			for className in self.states
			])


	def createProcesses(self):
		return collections.OrderedDict([
			(className, PROCESSES[className]())
			for className in self.processes
			])


	def createListeners(self):
		return collections.OrderedDict([
			(className, LISTENERS[className]())
			for className in self.listeners
			])


	def createLoggers(self):
		loggers = collections.OrderedDict()

		if self.logToShell:
			loggers["Shell"] = wholecell.loggers.shell.Shell(
				self.shellColumnHeaders
				)

		if self.logToDisk:
			loggers["Disk"] = wholecell.loggers.disk.Disk(
				self.outputDir,
				self.overwriteExistingFiles,
				self.logToDiskEvery
				)

		return loggers


	def createHooks(self):
		return collections.OrderedDict([
			(className, HOOKS[className]())
			for className in self.hooks
			])


	def toDict(self):
		return {
			attribute:getattr(self, attribute)
			for attribute in SIM_KWARG_DEFAULTS
			}

	# TODO: save, load

def getSimOptsFromEnvVars(optionsToNotGetFromEnvVars = None):
	optionsToNotGetFromEnvVars = optionsToNotGetFromEnvVars or []

	# We use this to check if any undefined WC_* environmental variables
	# were accidentally specified by the user
	wcEnvVars = [x for x in os.environ if x.startswith("WC_")]

	optionsAndEnvVars = dict(
		seed = ("WC_SEED", int),
		states = ("WC_STATES", json.loads),
		processes = ("WC_PROCESSES", json.loads),
		listeners = ("WC_LISTENERS", json.loads),
		hooks = ("WC_HOOKS", json.loads),
		lengthSec = ("WC_LENGTHSEC", int),
		timeStepSec = ("WC_TIMESTEPSEC", float),
		logToShell = ("WC_LOGTOSHELL", json.loads),
		shellColumnHeaders = ("WC_SHELLCOLUMNSHEADERS", json.loads),
		logToDisk = ("WC_LOGTODISK", json.loads),
		outputDir = ("WC_OUTPUTDIR", json.loads),
		overwriteExistingFiles = ("WC_OVERWRITEEXISTINGFILES", json.loads),
		logToDiskEvery = ("WC_LOGTODISKEVERY", int),
		kbLocation = ("WC_KBLOCATION", json.loads)
		)

	# These are options that the calling routine might set itself
	# While it could just overwrite them silently, removing them here
	# will at least alert the user
	for opt in optionsToNotGetFromEnvVars:
		del optionsAndEnvVars[opt]

	simOpts = {}

	# Get simulation options from environmental variables
	for option, (envVar, handler) in optionsAndEnvVars.iteritems():
		if os.environ.has_key(envVar) and len(os.environ[envVar]):
			simOpts[option] = handler(os.environ[envVar])
			wcEnvVars.remove(envVar)
		else:
			simOpts[option] = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS[option]
			if os.environ.has_key(envVar) and len(os.environ[envVar]) == 0:
				wcEnvVars.remove(envVar)

	# Check for extraneous environmental variables (probably typos by the user)
	assert (len(wcEnvVars) == 0), (
		"The following WC_* environmental variables were specified but " +
		"have no defined function: %s" % wcEnvVars
		)

	return simOpts