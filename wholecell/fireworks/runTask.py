"""
A command line to invoke a WCM Firetask with the given argument dict.
This is used to run a Firetask in a Docker container.
"""

from __future__ import absolute_import, division, print_function

import sys
import json

from wholecell.fireworks.firetasks import (
	ParcaTask,
	VariantSimDataTask,
	SimulationTask,
	SimulationDaughterTask,
	AnalysisVariantTask,
	AnalysisCohortTask,
	AnalysisSingleTask,
	AnalysisMultiGenTask,
	BuildCausalityNetworkTask,
	WriteJsonTask,
	)


FIRETASKS = (
	ParcaTask,
	VariantSimDataTask,
	SimulationTask,
	SimulationDaughterTask,
	AnalysisVariantTask,
	AnalysisCohortTask,
	AnalysisSingleTask,
	AnalysisMultiGenTask,
	BuildCausalityNetworkTask,
	WriteJsonTask,
	)
TASKS = {task.__name__: task for task in FIRETASKS}


if __name__ == '__main__':
	# Run the named WCM Firetask with the argument dict provided in JSON format.
	task_name = sys.argv[1]
	args = json.loads(sys.argv[2])
	assert isinstance(args, dict)

	print('runTask {}({})'.format(task_name, str(args)[1:-1]))

	task = TASKS[task_name](**args)
	task.run_task({})
