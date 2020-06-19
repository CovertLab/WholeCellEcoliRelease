from __future__ import absolute_import, division, print_function

from subprocess import call

from vivarium.core.composition import (
	make_agents,
	make_experiment_from_compartment_dicts,
	simulate_experiment,
	save_timeseries,
	plot_agents_multigen,
)
from vivarium.compartments.lattice import Lattice
from vivarium.plots.multibody_physics import plot_snapshots

from wholecell.utils.make_media import Media
from wholecell.states.external_state import ExternalState
from wcecoli_process import wcEcoliAgent
from wcecoli_compartment import WcEcoliCell


MEDIA_IDS = ["minimal", "minimal_minus_oxygen",
		"minimal_plus_amino_acids"]
MEDIA_ID = "minimal"
N_WCECOLI_AGENTS = 2
BOUNDS = (50, 50)
N_BINS = (20, 20)
OUT_DIR = 'out/agent'



def main():
	agent = wcEcoliAgent({})
	external_states = agent.ecoli_simulation.external_states
	# Assert agent has media_id MEDIA_ID
	recipe = external_states['Environment'].saved_media[MEDIA_ID]

	compartment = WcEcoliCell()
	agent_ids = ['wcecoli_{}'.format(i) for i in range(N_WCECOLI_AGENTS)]
	agents_dict = make_agents(agent_ids, compartment)

	environment_config = {
		'multibody': {
			'bounds': BOUNDS,
			'size': BOUNDS,
			'agents': agents_dict,
		},
		'diffusion': {
			'bounds': BOUNDS,
			'n_bins': N_BINS,
			'molecules': recipe.keys(),
            'depth': 10000.0,  # Deep to avoid depleting local molecules
			'diffusion': 5,  # 10x faster than the default 5e-1
			'gradient': {
				'type': 'linear',
				'molecules': {
					molecule: {
						'center': (0, 0),
						'base': concentration,
						'slope': 0,
					}
					for molecule, concentration in recipe.items()
				},
			},
		},
	}
	environment = Lattice(environment_config)
	environment_dict = environment.generate()

	emitter_dict = {'type': 'timeseries'}
	initial_state = {}
	experiment = make_experiment_from_compartment_dicts(
		environment_dict, agents_dict, emitter_dict, initial_state)
	settings = {
		'timestep': 1.0,
		'total_time': 20,
		'return_raw_data': True,
	}
	data = simulate_experiment(experiment, settings)

	# Plot snapshots
	agents = {
		time: time_data['agents'] for time, time_data in data.items()}
	fields = {
		time: time_data['fields'] for time, time_data in data.items()}
	snapshots_data = {
		'agents': agents,
		'fields': fields,
		'config': environment_config['multibody'],
	}
	plot_config = {
		'out_dir': OUT_DIR,
		'filename': 'snapshot',
	}
	plot_snapshots(snapshots_data, plot_config)

	# Plot Emitted Values
	plot_settings = {
		'agents_key': 'agents',
	}
	plot_agents_multigen(data, plot_settings, OUT_DIR, 'simulation')


if __name__ == '__main__':
	main()
