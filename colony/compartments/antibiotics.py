from __future__ import absolute_import, division, print_function

from vivarium.core.process import Generator
from vivarium.library.units import units
from vivarium.processes.derive_concentrations import (
	DeriveConcentrations)
from vivarium.processes.timeline import TimelineProcess
from vivarium_cell.processes.antibiotic_transport import AntibioticTransport
from vivarium_cell.processes.death import DeathFreezeState
from vivarium_cell.processes.derive_globals import DeriveGlobals
from vivarium_cell.processes.diffusion_cell_environment_ficks import (
	CellEnvironmentDiffusionFicks,
)

from colony.processes.antibiotic_hydrolysis import AntibioticHydrolysis
from colony.processes.wcecoli import WcEcoli
from colony.processes.wcecoli_derive_shape import WcEcoliDeriveShape
from colony.processes.wcecoli_meta_division import WcEcoliMetaDivision


INITIAL_INTERNAL_ANTIBIOTIC = 0
INITIAL_EXTERNAL_ANTIBIOTIC = 0
ANTIBIOTIC_KEY = 'nitrocefin'
PUMP_KEY = 'TRANS-CPLX-201[s]'
PORIN_KEY = 'porin'
# Source: (Wülfing & Plückthun, 1994)
PERIPLASM_FRACTION = 0.3
BETA_LACTAMASE_KEY = 'EG10040-MONOMER[p]'


class AntibioticsCell(Generator):
	'''Integrate antibiotic resistance and susceptibility with wcEcoli

	Integrates the WcEcoli process, which wraps the wcEcoli model, with
	processes to model antibiotic susceptibility (diffusion-based
	import and death) and resistance (hydrolysis and transport-based
	efflux). Also includes derivers.
	'''

	defaults = {
		'boundary_path': ('boundary',),
		'agents_path': ('..', '..', 'agents'),
		'fields_path': ('..', '..', 'fields'),
		'dimensions_path': ('..', '..', 'dimensions'),
		'update_fields': True,
		'antibiotic_transport': {
			'initial_pump': 0.0,
			'initial_internal_antibiotic': INITIAL_INTERNAL_ANTIBIOTIC,
			'intial_external_antibiotic': INITIAL_EXTERNAL_ANTIBIOTIC,
			# Reported in (Nagano & Nikaido, 2009)
			'pump_kcat': 1e1,  # Units: 1/s
			# Reported in (Nagano & Nikaido, 2009)
			'pump_km': 4.95e-3,  # Units: mM
			'pump_key': PUMP_KEY,
			'antibiotic_key': ANTIBIOTIC_KEY,
		},
		'hydrolysis': {
			'initial_catalyst': 0.0,
			'catalyst': BETA_LACTAMASE_KEY,
			'initial_target_internal': INITIAL_INTERNAL_ANTIBIOTIC,
			'target': ANTIBIOTIC_KEY,
			# Reported in (Galleni et al., 1988)
			'kcat': 490 / units.sec,
			# Reported in (Galleni et al., 1988)
			'Km': 500 * units.micromolar,
		},
		'death': {
			'detectors': {
				'antibiotic': {
					# (Kojima & Nikaido, 2013) reports values of 0.36
					# and 1.7 micro-molar for other beta-lactams
					# (benzylpenicillin and ampicillin), so we use 1 as
					# an estimate. Also note that (Nagano & Nikaido,
					# 2009) reports a value of 0.5 micro-molar for
					# nitrocefin from a reference, but we were unable to
					# find where the reference reports the value.
					'antibiotic_threshold': 1e-3,  # Units: mM
					'antibiotic_key': ANTIBIOTIC_KEY,
				},
			},
			'targets': [
				'wcecoli', 'meta_division', 'death',
				'antibiotic_transport'
			],
		},
		'derive_shape': {},
		'cell_environment_diffusion': {
			'default_state': {
				'external': {
					ANTIBIOTIC_KEY: INITIAL_EXTERNAL_ANTIBIOTIC,
				},
				'internal': {
					ANTIBIOTIC_KEY: INITIAL_INTERNAL_ANTIBIOTIC,
				},
				'global': {
					'periplasm_volume': (
						1.2 * units.fL * PERIPLASM_FRACTION),
				},
			},
			'molecules_to_diffuse': [ANTIBIOTIC_KEY],
			# (Nagano & Nikaido, 2009) reports that their mutant strain,
			# RAM121, has 10-fold faster influx of nitrocefin with a
			# permeability of 0.2e-5 cm/s, so wildtype has a
			# permeability of 0.2e-6 cm/s.
			'permeability': 0.2e-6 * units.cm / units.sec,
			# From (Nagano & Nikaido, 2009)
			'surface_area_mass_ratio': 132 * units.cm**2 / units.mg,
			'volume_variable': 'periplasm_volume',
		},
		'derive_concentrations': {
			'concentration_keys': [PUMP_KEY, BETA_LACTAMASE_KEY],
		},
		'timeline': {},
		'derive_globals': {
			'periplasm_volume_fraction': PERIPLASM_FRACTION,
		},
	}

	def generate_processes(self, config):
		wcecoli_process = WcEcoli(config)
		meta_division = WcEcoliMetaDivision({
			'agent_id': config['agent_id'],
			'compartment': self,
		})
		derive_shape = WcEcoliDeriveShape(config['derive_shape'])
		death = DeathFreezeState(config['death'])
		antibiotic_transport = AntibioticTransport(
			config['antibiotic_transport'])
		hydrolysis = AntibioticHydrolysis(
			config['hydrolysis'])
		cell_environment_diffusion = CellEnvironmentDiffusionFicks(
			config['cell_environment_diffusion'])
		derive_concentrations = DeriveConcentrations(
			config['derive_concentrations'])
		timeline = TimelineProcess(config['timeline'])
		derive_globals = DeriveGlobals(config['derive_globals'])
		return {
			'derive_globals': derive_globals,
			'derive_concentrations': derive_concentrations,
			'wcecoli': wcecoli_process,
			'meta_division': meta_division,
			'derive_shape': derive_shape,
			'death': death,
			'antibiotic_transport': antibiotic_transport,
			'hydrolysis': hydrolysis,
			'cell_environment_diffusion': cell_environment_diffusion,
			'timeline': timeline,
		}

	def generate_topology(self, config=None):
		boundary_path = config['boundary_path']
		topology = {
			'derive_globals': {
				'global': boundary_path,
			},
			'derive_concentrations': {
				'global': boundary_path,
				'counts': boundary_path + ('bulk_molecules_report',),
				'concentrations': (
					boundary_path + ('bulk_molecule_concentrations',)
				),
			},
			'wcecoli': {
				'bulk_molecules_report': (
					boundary_path + ('bulk_molecules_report',)
				),
				'unique_molecules_report': (
					boundary_path + ('unique_molecules_report',)
				),
				'listeners_report': (
					boundary_path + ('listeners_report',)
				),
				'global': boundary_path,
				'external': (
					boundary_path + ('external',)
				),
				'dimensions': config['dimensions_path'],
			},
			'meta_division': {
				'global': boundary_path,
				'cells': config['agents_path'],
			},
			'derive_shape': {
				'global': boundary_path,
			},
			'death': {
				'global': boundary_path,
				'internal': boundary_path + ('cytoplasm',),
				'processes': tuple(),
			},
			'antibiotic_transport': {
				'internal': boundary_path + ('cytoplasm',),
				'external': boundary_path + ('external',),
				'pump_port': (
					boundary_path + ('bulk_molecule_concentrations',)
				),
				'fields': config['fields_path'],
				'fluxes': boundary_path + ('fluxes',),
				'global': boundary_path,
				'dimensions': config['dimensions_path'],
			},
			'hydrolysis': {
				'internal': boundary_path + ('cytoplasm',),
				'catalyst_port': (
					boundary_path + ('bulk_molecule_concentrations',)
				),
				'fields': config['fields_path'],
				'fluxes': boundary_path + ('fluxes',),
				'global': boundary_path,
				'dimensions': config['dimensions_path'],
			},
			'cell_environment_diffusion': {
				'internal': boundary_path + ('cytoplasm',),
				'external': boundary_path + ('external',),
				'fields': config['fields_path'],
				'global': boundary_path,
				'dimensions': config['dimensions_path'],
			},
			'timeline': {
				'global': boundary_path,
			},
		}
		if len(config['timeline']['timeline']) > 1:
			# Then we change the fields
			topology['timeline']['fields'] = config['fields_path']
		if config['update_fields']:
			topology['wcecoli']['fields'] = config['fields_path']
		return topology
