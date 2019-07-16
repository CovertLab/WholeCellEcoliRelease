#!/usr/bin/env python

"""
Create a demo workflow for quick tests of the workflow software.
"""

from __future__ import absolute_import, division, print_function

import argparse
import os

import wholecell.utils.filepath as fp
from runscripts.sisyphus.util.workflow import Task, Workflow


DOCKER_IMAGE = 'python:2.7.16'
STORAGE_PREFIX_ROOT = 'sisyphus:data/'


def demo(worker_count=0):
	owner_id = os.environ['USER']
	timestamp = fp.timestamp()
	name = 'Demo_{}_{}'.format(owner_id, timestamp)
	storage_prefix = os.path.join(STORAGE_PREFIX_ROOT, owner_id, 'demo', timestamp, '')
	wf = Workflow(name)

	filename = '/tmp/output.log'
	code = ("with open('" + filename + "', 'w') as f:\n"
        "  for i in range(100):\n"
        "    f.write('This is line {}\\n'.format(i))\n"
        "    print 'hello', i")
	line_task = Task((),
		name='lines',
		image=DOCKER_IMAGE,
		inputs=[],
		outputs=[filename],
		storage_prefix=storage_prefix,
		local_prefix='/tmp',
		commands=[{'command': ['python', '-u', '-c', code]}])
	wf.add_task(line_task)

	wf.send(worker_count=worker_count)

def cli():
	parser = argparse.ArgumentParser(description='Demo workflow')
	parser.add_argument('-w', '--workers', type=int, default=0,
		help='number of worker nodes to launch; default = 0')
	args = parser.parse_args()
	demo(worker_count=args.workers)


if __name__ == '__main__':
	cli()
