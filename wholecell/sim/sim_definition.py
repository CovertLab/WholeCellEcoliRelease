"""

Defines a class for storing and managing the parameters involved in simulation
instantiation.

"""

from __future__ import division

import collections

# import tables

# References to sub-simulation abstractions

# States
import wholecell.states.bulk_molecules
import wholecell.states.unique_molecules
import wholecell.states.chromosome
import wholecell.states.transcripts

STATE_CLASSES = (
	wholecell.states.bulk_molecules.BulkMolecules,
	wholecell.states.unique_molecules.UniqueMolecules,
	wholecell.states.chromosome.Chromosome,
	wholecell.states.transcripts.Transcripts,
	)

STATES = {stateClass.name():stateClass for stateClass in STATE_CLASSES}

# Processes
import wholecell.processes.complexation
import wholecell.processes.metabolism
import wholecell.processes.metabolism_fba
import wholecell.processes.rna_degradation
import wholecell.processes.transcription.bulk_transcription
import wholecell.processes.translation.translation
import wholecell.processes.transcription.toy_transcription
import wholecell.processes.replication
import wholecell.processes.transcription.transcription_net
import wholecell.processes.translation.translation_net
import wholecell.processes.translation.unique_polypeptide_initiation
import wholecell.processes.translation.unique_polypeptide_elongation
import wholecell.processes.translation.unique_polypeptide_elongation_enzlim
import wholecell.processes.transcription.unique_transcript_initiation
import wholecell.processes.transcription.unique_transcript_elongation
import wholecell.processes.transcription.unique_transcript_elongation_enzlim
import wholecell.processes.protein_degradation

PROCESS_CLASSES = (
	wholecell.processes.complexation.Complexation,
	wholecell.processes.metabolism.Metabolism,
	wholecell.processes.metabolism_fba.MetabolismFba,
	wholecell.processes.rna_degradation.RnaDegradation,
	wholecell.processes.transcription.bulk_transcription.BulkTranscription,
	wholecell.processes.translation.translation.Translation,
	wholecell.processes.transcription.toy_transcription.ToyTranscription,
	wholecell.processes.replication.Replication,
	wholecell.processes.transcription.transcription_net.TranscriptionNet,
	wholecell.processes.translation.translation_net.TranslationNet,
	wholecell.processes.translation.unique_polypeptide_initiation.UniquePolypeptideInitiation,
	wholecell.processes.translation.unique_polypeptide_elongation.UniquePolypeptideElongation,
	wholecell.processes.translation.unique_polypeptide_elongation_enzlim.UniquePolypeptideElongationEnzlim,
	wholecell.processes.transcription.unique_transcript_initiation.UniqueTranscriptInitiation,
	wholecell.processes.transcription.unique_transcript_elongation.UniqueTranscriptElongation,
	wholecell.processes.transcription.unique_transcript_elongation_enzlim.UniqueTranscriptElongationEnzlim,
	wholecell.processes.protein_degradation.ProteinDegradation
	)

PROCESSES = {processClass.name():processClass for processClass in PROCESS_CLASSES}

# Listeners
import wholecell.listeners.mass
import wholecell.listeners.metabolic_flux
import wholecell.listeners.replication_fork
import wholecell.listeners.ntp_usage
import wholecell.listeners.aa_usage

LISTENER_CLASSES = (
	wholecell.listeners.mass.Mass,
	wholecell.listeners.metabolic_flux.MetabolicFlux,
	wholecell.listeners.replication_fork.ReplicationForkPosition,
	wholecell.listeners.ntp_usage.NtpUsage,
	wholecell.listeners.aa_usage.AAUsage
	)

LISTENERS = {listenerClass.name():listenerClass for listenerClass in LISTENER_CLASSES}

# Loggers
import wholecell.loggers.shell
import wholecell.loggers.disk

# TODO: logger logic more consistent with listeners/states/processes

# Hooks
# TODO: move hooks to their own directory/files
import wholecell.sim.hooks

HOOK_CLASSES = (
	wholecell.sim.hooks.RnapCountHook,
	wholecell.sim.hooks.RibosomeCountHook
	)

HOOKS = {hookClass.name():hookClass for hookClass in HOOK_CLASSES}

# Default parameters

DEFAULT_STATES = (
	'BulkMolecules',
	'UniqueMolecules'
	)

DEFAULT_PROCESSES = (
	'Metabolism',
	'RnaDegradation',
	'UniqueTranscriptInitiation',
	'UniqueTranscriptElongation',
	'UniquePolypeptideInitiation',
	'UniquePolypeptideElongation',
	'Replication',
	'ProteinDegradation'
	)

DEFAULT_LISTENERS = (
	'Mass',
	'ReplicationForkPosition',
	'NtpUsage',
	'AAUsage'
	)

DEFAULT_HOOKS = ( # NOTE: there should probably never be any default hooks
	)

DEFAULT_LENGTH = 3600 # sec
DEFAULT_TIME_STEP = 1 # sec

DEFAULT_SEED = None

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
	logToShell = True,
	logToDisk = False, outputDir = None, overwriteExistingFiles = False, logToDiskEvery = None,
	rebuildKB = True
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
		loggers = []

		if self.logToShell:
			loggers.append(
				wholecell.loggers.shell.Shell()
				)

		if self.logToDisk:
			loggers.append(
				wholecell.loggers.disk.Disk(
					self.outputDir,
					self.overwriteExistingFiles,
					self.logToDiskEvery
					)
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
