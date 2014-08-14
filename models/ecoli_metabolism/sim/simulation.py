
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

# Initialization
from models.ecoli_metabolism.sim.initial_conditions import calcInitialConditions


class EcoliMetabolismSimulation(Simulation):
	_stateClasses = (
		BulkMolecules,
		)

	_processClasses = (
		TranscriptElongation,
		# RnaDegradation,
		PolypeptideElongation,
		# ProteinDegradation,
		# Replication,
		Metabolism,
		# Maintenance,
		)

	_listenerClasses = (
		Mass,
		# NtpUsage, # restore these in some general sense...
		# AAUsage,
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
