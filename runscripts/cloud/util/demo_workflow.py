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
	def add_task(name='', inputs=(), outputs=(), command=()):
		task = Task(
			name=name,
			image=DOCKER_IMAGE,
			inputs=inputs,
			outputs=outputs,
			storage_prefix=storage_prefix,
			internal_prefix='/tmp',
			command=command)
		wf.add_task(task)

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
	add_task(
		name='lines',
		outputs=[lines_filename],
		command=['python', '-u', '-c', code])

	add_task(
		name='count',
		inputs=[lines_filename],
		outputs=['>/tmp/count.log'],
		command=['wc', lines_filename])

	add_task(
		name='error_test',
		inputs=[],  # test error handling by "forgetting" to download the input
		outputs=['>/tmp/error-test.log'],
		command=['wc', lines_filename])

	add_task(
		name='index_exception',
		outputs=['>/tmp/index_exception.log'],
		command=['python', '-u', '-c', "()[1]"])

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
