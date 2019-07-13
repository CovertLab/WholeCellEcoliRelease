from __future__ import absolute_import, division, print_function

import sys
import json

from wholecell.fireworks.firetasks import ParcaTask
from wholecell.fireworks.firetasks import VariantSimDataTask
from wholecell.fireworks.firetasks import SimulationTask
from wholecell.fireworks.firetasks import SimulationDaughterTask
from wholecell.fireworks.firetasks import AnalysisVariantTask
from wholecell.fireworks.firetasks import AnalysisCohortTask
from wholecell.fireworks.firetasks import AnalysisSingleTask
from wholecell.fireworks.firetasks import AnalysisMultiGenTask
from wholecell.fireworks.firetasks import BuildCausalityNetworkTask
from wholecell.fireworks.firetasks import WriteJsonTask


TASKS = {
	'parca': ParcaTask,
	'variant_sim_data': VariantSimDataTask,
	'simulation': SimulationTask,
	'simulation_daughter': SimulationDaughterTask,
	'analysis_variant': AnalysisVariantTask,
	'analysis_cohort': AnalysisCohortTask,
	'analysis_single': AnalysisSingleTask,
	'analysis_multigen': AnalysisMultiGenTask,
	'build_causality_network': BuildCausalityNetworkTask,
	'write_json': WriteJsonTask}

if __name__ == '__main__':
	# Run the named WCM Firetask with the argument dict provided in JSON format.
	task_name = sys.argv[1]
	args = json.loads(sys.argv[2])
	assert isinstance(args, dict)

	print('runTask {}({})'.format(task_name, str(args)[1:-1]))

	task = TASKS[task_name](**args)
	task.run_task({})
