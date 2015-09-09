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
from models.ecoli.processes.replication_elongation import ReplicationElongation
from models.ecoli.processes.polypeptide_initiation import PolypeptideInitiation
from models.ecoli.processes.polypeptide_elongation import PolypeptideElongation
from models.ecoli.processes.transcript_initiation import TranscriptInitiation
from models.ecoli.processes.transcript_elongation import TranscriptElongation
from models.ecoli.processes.protein_degradation import ProteinDegradation
from models.ecoli.processes.atp_usage import AtpUsage
from models.ecoli.processes.chromosome_formation import ChromosomeFormation
from models.ecoli.processes.equilibrium import Equilibrium

# Listeners
from models.ecoli.listeners.mass import Mass
from models.ecoli.listeners.replication_data import ReplicationData
from models.ecoli.listeners.ribosome_data import RibosomeData
from models.ecoli.listeners.gene_copy_number import GeneCopyNumber
from models.ecoli.listeners.unique_molecule_counts import UniqueMoleculeCounts
from models.ecoli.listeners.fba_results import FBAResults
from models.ecoli.listeners.rna_degradation_listener import RnaDegradationListener
from models.ecoli.listeners.transcript_elongation_listener import TranscriptElongationListener
from models.ecoli.listeners.rnap_data import RnapData
from models.ecoli.listeners.enzyme_kinetics import EnzymeKinetics
from models.ecoli.listeners.growth_limits import GrowthLimits


# Analysis
import models.ecoli.analysis.single
import models.ecoli.analysis.cohort

from models.ecoli.sim.initial_conditions import calcInitialConditions
from wholecell.sim.divide_cell import divide_cell
from models.ecoli.sim.initial_conditions import setDaughterInitialConditions

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
		ReplicationElongation,
		ProteinDegradation,
		Complexation,
		AtpUsage,
		ChromosomeFormation,
		Equilibrium,
		)

	_listenerClasses = (
		Mass,
		ReplicationData,
		RibosomeData,
		GeneCopyNumber,
		UniqueMoleculeCounts,
		FBAResults,
		RnaDegradationListener,
		TranscriptElongationListener,
		RnapData,
		EnzymeKinetics,
		GrowthLimits
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
		"Expected fold change"
		]

	_logToDisk = False

class EcoliDaughterSimulation(EcoliSimulation):
	_initialConditionsFunction = setDaughterInitialConditions
