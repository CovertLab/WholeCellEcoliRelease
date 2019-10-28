
from wholecell.sim.simulation import Simulation

# States
from wholecell.states.bulk_molecules import BulkMolecules
from wholecell.states.unique_molecules import UniqueMolecules
from wholecell.states.local_environment import LocalEnvironment

# Processes
from models.ecoli.processes.complexation import Complexation
from models.ecoli.processes.metabolism import Metabolism
from models.ecoli.processes.rna_degradation import RnaDegradation
from models.ecoli.processes.chromosome_replication import ChromosomeReplication
from models.ecoli.processes.polypeptide_initiation import PolypeptideInitiation
from models.ecoli.processes.polypeptide_elongation import PolypeptideElongation
from models.ecoli.processes.transcript_initiation import TranscriptInitiation
from models.ecoli.processes.transcript_elongation import TranscriptElongation
from models.ecoli.processes.protein_degradation import ProteinDegradation
from models.ecoli.processes.equilibrium import Equilibrium
from models.ecoli.processes.tf_binding import TfBinding
from models.ecoli.processes.two_component_system import TwoComponentSystem

# Listeners
from models.ecoli.listeners.mass import Mass
from models.ecoli.listeners.replication_data import ReplicationData
from models.ecoli.listeners.ribosome_data import RibosomeData
from models.ecoli.listeners.unique_molecule_counts import UniqueMoleculeCounts
from models.ecoli.listeners.fba_results import FBAResults
from models.ecoli.listeners.rna_degradation_listener import RnaDegradationListener
from models.ecoli.listeners.transcript_elongation_listener import TranscriptElongationListener
from models.ecoli.listeners.rnap_data import RnapData
from models.ecoli.listeners.enzyme_kinetics import EnzymeKinetics
from models.ecoli.listeners.growth_limits import GrowthLimits
from models.ecoli.listeners.cell_division import CellDivision
from models.ecoli.listeners.rna_synth_prob import RnaSynthProb
from models.ecoli.listeners.monomer_counts import MonomerCounts
from models.ecoli.listeners.mRNA_counts import mRNACounts
from models.ecoli.listeners.complexation_listener import ComplexationListener
from models.ecoli.listeners.equilibrium_listener import EquilibriumListener

from models.ecoli.sim.initial_conditions import calcInitialConditions
from wholecell.sim.divide_cell import divide_cell
from models.ecoli.sim.initial_conditions import setDaughterInitialConditions

class EcoliSimulation(Simulation):
	_internalStateClasses = (
		BulkMolecules,
		UniqueMolecules,
		)

	_externalStateClasses = (
		LocalEnvironment,
		)

	_first_process_classes = (
		RnaDegradation,
		TranscriptInitiation,
		TranscriptElongation,
		PolypeptideInitiation,
		PolypeptideElongation,
		ChromosomeReplication,
		ProteinDegradation,
		Complexation,
		Equilibrium,
		TfBinding,
		TwoComponentSystem,
		)

	_second_process_classes = (
		Metabolism,
		)

	_processClasses = _first_process_classes + _second_process_classes

	_listenerClasses = (
		Mass,
		ReplicationData,
		RibosomeData,
		UniqueMoleculeCounts,
		FBAResults,
		RnaDegradationListener,
		TranscriptElongationListener,
		RnapData,
		EnzymeKinetics,
		GrowthLimits,
		CellDivision,
		RnaSynthProb,
		MonomerCounts,
		mRNACounts,
		ComplexationListener,
		EquilibriumListener,
		)

	_hookClasses = ()

	_initialConditionsFunction = calcInitialConditions

	_divideCellFunction = divide_cell

	_logToShell = True
	_shellColumnHeaders = [
		"Time (s)",
		"Dry mass (fg)",
		"Dry mass fold change",
		"Protein fold change",
		"RNA fold change",
		"Small mol fold change",
		"Expected fold change"
		]

	_logToDisk = False

class EcoliDaughterSimulation(EcoliSimulation):
	_initialConditionsFunction = setDaughterInitialConditions


def ecoli_simulation(**options):
	"""Instantiate an initial EcoliSimulation or a daughter
	EcoliDaughterSimulation with the given options, depending on whether
	there's a non-None `inheritedStatePath` option.
	"""
	is_daughter = options.get('inheritedStatePath', None) is not None
	return EcoliDaughterSimulation(**options) if is_daughter else EcoliSimulation(**options)
