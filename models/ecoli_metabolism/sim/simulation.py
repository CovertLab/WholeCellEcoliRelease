
import os

from wholecell.sim.simulation import Simulation

# States
from wholecell.states.bulk_molecules import BulkMolecules

# Processes

# Listeners
from models.ecoli.listeners.mass import Mass


class EcoliMetabolismSimulation(Simulation):
	_stateClasses = (
		BulkMolecules,
		)

	_processClasses = (
		# TranscriptInitiation,
		# TranscriptElongation,
		# RnaDegradation,
		# PolypeptideInitiation,
		# PolypeptideElongation,
		# ProteinDegradation,
		# Replication,
		# Metabolism,
		# AtpUsage
		)

	_listenerClasses = (
		Mass,
		# NtpUsage, # restore these in some general sense...
		# AAUsage,
		)

	_hookClasses = ()

	_initialConditionsFunction = lambda sim, kb: None

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
