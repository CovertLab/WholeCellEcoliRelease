"""

Defines a class for storing and managing the parameters involved in simulation
instantiation.

"""

from __future__ import division

import collections

# import tables

import wholecell.states.mass
import wholecell.states.bulk_molecules
import wholecell.states.unique_molecules
import wholecell.states.chromosome
import wholecell.states.transcripts

STATE_CLASSES = [
	wholecell.states.mass.Mass,
	wholecell.states.bulk_molecules.BulkMolecules,
	wholecell.states.unique_molecules.UniqueMolecules,
	wholecell.states.chromosome.Chromosome,
	wholecell.states.transcripts.Transcripts,
	]

import wholecell.processes.complexation
import wholecell.processes.metabolism
import wholecell.processes.metabolism_fba
import wholecell.processes.rna_degradation
import wholecell.processes.transcription.bulk_transcription
import wholecell.processes.translation.translation
import wholecell.processes.free_production
import wholecell.processes.transcription.toy_transcription
import wholecell.processes.toy_protein_degradation
import wholecell.processes.toy_replication
import wholecell.processes.replication
import wholecell.processes.transcription.transcription_net
import wholecell.processes.translation.translation_net
import wholecell.processes.translation.unique_polypeptide_initiation
import wholecell.processes.translation.unique_polypeptide_elongation
import wholecell.processes.transcription.unique_transcript_initiation
import wholecell.processes.transcription.unique_transcript_elongation

PROCESS_CLASSES = [
	wholecell.processes.metabolism.Metabolism,
	wholecell.processes.metabolism_fba.MetabolismFba,
	wholecell.processes.rna_degradation.RnaDegradation,
	wholecell.processes.transcription.bulk_transcription.BulkTranscription,
	wholecell.processes.translation.translation.Translation,
	wholecell.processes.free_production.FreeProduction,
	wholecell.processes.transcription.toy_transcription.ToyTranscription,
	wholecell.processes.toy_protein_degradation.ToyProteinDegradation,
	wholecell.processes.toy_replication.ToyReplication,
	wholecell.processes.replication.Replication,
	wholecell.processes.transcription.transcription_net.TranscriptionNet,
	wholecell.processes.translation.translation_net.TranslationNet,
	wholecell.processes.translation.unique_polypeptide_initiation.UniquePolypeptideInitiation,
	wholecell.processes.translation.unique_polypeptide_elongation.UniquePolypeptideElongation,
	wholecell.processes.transcription.unique_transcript_initiation.UniqueTranscriptInitiation,
	wholecell.processes.transcription.unique_transcript_elongation.UniqueTranscriptElongation
	]

STATES = {stateClass.name():stateClass for stateClass in STATE_CLASSES}
PROCESSES = {processClass.name():processClass for processClass in PROCESS_CLASSES}

DEFAULT_STATES = [
	'Mass',
	'BulkMolecules',
	'UniqueMolecules'
	]

DEFAULT_PROCESSES = [
	'Metabolism',
	'RnaDegradation',
	'UniqueTranscriptInitiation',
	'UniqueTranscriptElongation',
	'UniquePolypeptideInitiation',
	'UniquePolypeptideElongation',
	'Replication'
	]

SIM_KWARG_DEFAULTS = dict(
	includedStates = DEFAULT_STATES, includedProcesses = DEFAULT_PROCESSES,
	lengthSec = None, timeStepSec = None,
	seed = None,
	reconstructKB = False,
	logToShell = True,
	logToDisk = False, outputDir = None, overwriteExistingFiles = False, logToDiskEvery = None
	)

class SimDefinitionException(Exception): pass

# TODO: incorporate/define more of the default simulation behavior here
# need to determine what parts of simulation
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
			(className, STATE_CLASSES[className])()
			for className in self.includedStates
			])


	def createProcesses(self):
		return collections.OrderedDict([
			(className, PROCESS_CLASSES[className])()
			for className in self.includedProcesses
			])

	# TODO: save, load
