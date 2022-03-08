#! /usr/bin/env python
"""
Script to rerun fireworks that get preempted when submitted to the owners node
on sherlock.  Preempted jobs do not have a way of updating the status in mongoDB
to FIZZLED and will be stuck on RUNNING.  This checks the name of RUNNING tasks
against those in slurm to see if any RUNNING tasks are not actually running.

Use the option -f to automatically set these tasks to rerun, otherwise only the
fw_ids will be printed for manually rerunning them.

TODO:
- add a way of catching sim issues that cause tasks to be stuck on RUNNING
(MAX_RERUN is used to prevent the same task from being rerun many times if the
issue comes from the code)
- parse all options from command line for more flexibility
"""

import json
import os
import sys
import time
from typing import Dict

from wholecell.utils import filepath


VERBOSE = False
MAX_RERUN = 5
CACHE_WAIT = 60
SLEEP = 300
RUN_CMD = len(sys.argv) > 1 and sys.argv[1] == '-f'


def get_fws(criteria):
	result = json.loads(filepath.run_cmdline(f'lpad get_fws {criteria}', fallback='[]'))

	# Handle single fw case to always make a list of dicts
	if isinstance(result, dict):
		result = [result]
	return result

rerun_counts = {}  # type: Dict[str, int]
while True:
	# Get all fireworks with running status and all active job names
	running = get_fws('-s RESERVED') + get_fws('-s RUNNING')
	if VERBOSE:
		print(running)
	time.sleep(CACHE_WAIT)  # squeue is sometimes cached so wait until we know the new tasks will be returned
	if (squeue := filepath.run_cmdline(f'squeue -u {os.environ["USER"]} -o %j')) is None:
		continue
	else:
		tasks = set(squeue.split('\n'))
	if VERBOSE:
		print(tasks)

	# Check if running fireworks are not in the slurm queue
	fw_ids = ','.join([str(fw['fw_id']) for fw in running if fw['name'] not in tasks])
	if VERBOSE:
		print(fw_ids)

	if fw_ids:
		# Get the current status of fireworks in case they finished between
		# running get_fws and squeue calls
		current = get_fws(f'-i {fw_ids}')
		if VERBOSE:
			print(current)

		# Rerun fireworks that still have a running status but are not in slurm
		fws = []
		for fw in current:
			if fw['state'] == 'RUNNING' or fw['state'] == 'RESERVED':
				fw_id = str(fw['fw_id'])
				rerun_counts[fw_id] = rerun_counts.get(fw_id, 0) + 1
				if rerun_counts[fw_id] <= MAX_RERUN:
					fws.append(fw_id)
		fws_to_rerun = ','.join(fws)
		if fws_to_rerun:
			if RUN_CMD:
				print(f'Rerunning fw_ids {fws_to_rerun}')
				result = filepath.run_cmdline(f'lpad rerun_fws -i {fws_to_rerun}', input_='y', timeout=None)
				if result is None:
					filepath.run_cmdline('lpad admin unlock')
			else:
				print(f'Would rerun fw_ids {fws_to_rerun} if -f passed as arg. Manually rerun with:'
					f'\n\tlpad rerun_fws -i {fws_to_rerun}')

	# Sleep for a bit
	print(f'{time.ctime()}: sleeping for {SLEEP} s...')
	time.sleep(SLEEP)
