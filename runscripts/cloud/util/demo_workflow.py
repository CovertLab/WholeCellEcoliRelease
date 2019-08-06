#!/usr/bin/env python

"""
Create a demo workflow for quick tests of the workflow software.
"""

from __future__ import absolute_import, division, print_function

import argparse
import os

import wholecell.utils.filepath as fp
from runscripts.cloud.util.workflow import Task, Workflow, STDOUT_PATH


DOCKER_IMAGE = 'python:2.7.16'
STORAGE_PREFIX_ROOT = 'sisyphus:data/'


def demo(worker_count=0, dump=False):
	owner_id = os.environ['USER']
	timestamp = fp.timestamp()
	name = 'Demo_{}_{}'.format(owner_id, timestamp)
	storage_prefix = os.path.join(STORAGE_PREFIX_ROOT, owner_id, 'demo', timestamp, '')
	wf = Workflow(name)

	# Build the workflow.
	lines_filename = '/tmp/lines.txt'
	code = ("with open('" + lines_filename + "', 'w') as f:\n"
        "  for i in range(100):\n"
        "    f.write('This is line {}\\n'.format(i))\n"
        "    print 'hello', i")
	line_task = Task(
		name='lines',
		image=DOCKER_IMAGE,
		outputs=[lines_filename],
		storage_prefix=storage_prefix,
		internal_prefix='/tmp',
		command=['python', '-u', '-c', code])
	wf.add_task(line_task)

	count_task = Task(
		name='count',
		image=DOCKER_IMAGE,
		inputs=[lines_filename],
		outputs=['>/tmp/count.txt'],
		storage_prefix=storage_prefix,
		internal_prefix='/tmp',
		command=['wc', lines_filename])
	wf.add_task(count_task)

	# Dump or run the workflow.
	if dump:
		wf.write()
	else:
		wf.send(worker_count=worker_count)

def cli():
	parser = argparse.ArgumentParser(description='Demo workflow')
	parser.add_argument('-w', '--workers', type=int, default=1,
		help='number of worker nodes to launch; default = 1')
	parser.add_argument('--dump', action='store_true',
		help='Dump the built workflow to JSON files for your review *instead* of'
			 ' sending them to the Gaia workflow server. This is useful for'
			 ' testing and debugging. You can upload them manually or re-run'
			 ' this program without `--dump`.')

	args = parser.parse_args()
	demo(worker_count=args.workers, dump=args.dump)


if __name__ == '__main__':
	cli()
