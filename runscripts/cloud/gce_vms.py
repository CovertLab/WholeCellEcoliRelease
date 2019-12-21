#!/usr/bin/env python
"""Launch Google Compute Engine VMs."""

from __future__ import absolute_import, division, print_function

import argparse
import os
import re
import sys

if os.name == 'posix' and sys.version_info[0] < 3:
	import subprocess32 as subprocess
else:
	import subprocess

from typing import Any, Dict, Iterable, List, Optional
from wholecell.utils import filepath as fp


MAX_VMS = 100  # don't launch more than this many GCE VMs at a time

SERVICE_ACCOUNT = '441871726775-compute@developer.gserviceaccount.com'

SCOPES = ('https://www.googleapis.com/auth/devstorage.read_only,'
		  'https://www.googleapis.com/auth/logging.write,'
		  'https://www.googleapis.com/auth/monitoring.write,'
		  'https://www.googleapis.com/auth/servicecontrol,'
		  'https://www.googleapis.com/auth/service.management.readonly,'
		  'https://www.googleapis.com/auth/trace.append')


def get_gloud_get(section_property):
	# type: (str) -> str
	"""Get a `section/property` configuration setting from the gcloud command line tool."""
	value = fp.run_cmd(['gcloud', 'config', 'get-value', str(section_property)])
	assert value, 'You must configure `gcloud config set {} VALUE`'.format(section_property)
	return value


def get_gloud_project():
	# type: () -> str
	"""Get the current Google Cloud Platform (GCP) project from gcloud."""
	return get_gloud_get('core/project')


def get_gloud_zone():
	# type: () -> str
	"""Get the current Google Compute Engine (GCE) zone from gcloud."""
	return get_gloud_get('compute/zone')


def make_VM_names(name_prefix, base=0, count=1):
	# type: (str, int, int) -> List[str]
	"""Make a list of valid GCE VM names from the prefix, base index, and count."""
	sanitized = re.sub(
		r'[^-a-z0-9]+', '-',
		name_prefix.lower().replace('workflow', ''))
	names = ['{}-{}'.format(sanitized, i) for i in range(base, base + count)]
	return names


def launch_GCE_VMs(name_prefix, base=0, count=1, command_options=None, **metadata):
	# type: (str, int, int, Optional[Dict[str, Any]], **Any) -> None
	"""Parallel-launch VMs on GCE.

	This provides default command options for `gcloud compute instances create`.
	The caller should at least set the `image-family` option and override
	defaults as needed.

	This converts command_options and metadata to strings, then passes tokens
	to `gcloud` avoiding quoting problems, but the metadata fields get joined
	with commas so those keys and values must not contain commas.
	"""
	assert name_prefix, 'the name_prefix must not be empty'
	assert 0 <= count < MAX_VMS, 'VM count ({}) must be in the range [0 .. {}]'.format(
		count, MAX_VMS)
	instance_names = make_VM_names(name_prefix, base, count)

	print('Launching {} Google Compute Engine VM(s).'.format(count))
	if count <= 0:
		return

	project = get_gloud_project()
	options = {
		'project': project,
		'image-project': project,
		'zone': get_gloud_zone(),
		'machine-type': 'n1-standard-1',
		'subnet': 'default',
		'network-tier': 'PREMIUM',
		'maintenance-policy': 'MIGRATE',
		'boot-disk-size': '200GB',
		'boot-disk-type': 'pd-standard',
		'service-account': SERVICE_ACCOUNT,
		'scopes': SCOPES,
		}
	options.update(command_options or {})

	metadata_string = ','.join('{}={}'.format(k, v) for k, v in metadata.items())
	if metadata_string:
		options['metadata'] = metadata_string

	options_list = ['--{}={}'.format(k, v) for k, v in options.items()]
	cmd_tokens = ['gcloud', 'compute', 'instances', 'create'
		] + instance_names + options_list
	subprocess.call(cmd_tokens, env=os.environ)


def launch_sisyphus_workers(name_prefix, workflow, base=0, count=1):
	# type: (str, str, int, int) -> None
	"""Parallel-launch Sisyphus worker VMs on GCE for the given workflow name."""
	assert workflow, 'a workflow name is required'

	options = {'image-family': 'sisyphus-worker'}
	launch_GCE_VMs(name_prefix, base, count, options, workflow=workflow,
		description='sisyphus worker')


def launch_fireworkers(name_prefix, base=0, count=1, **metadata):
	# type: (str, int, int, **Any) -> None
	"""Parallel-launch FireWorks worker VMs on GCE."""
	options = {'image-family': 'fireworker'}
	launch_GCE_VMs(name_prefix, base, count, options, description='fireworker',
		**metadata)


def main():
	parser = argparse.ArgumentParser(
		description='Launch Google Compute Engine VMs all at once.')
	parser.add_argument('-s', '--sisyphus', action='store_true',
		help='Launch Sisyphus workers (rather than Fireworkers)')
	parser.add_argument('-w', '--workflow',
		help='the workflow name; required with --sisyphus')
	parser.add_argument('name_prefix', metavar='NAME-PREFIX',
		help="the GCE VM name prefix; it'll get numeric suffixes")
	parser.add_argument('-i', '--index', type=int, default=0,
		help='the base index number for the VM names numeric suffixes')
	parser.add_argument('-c', '--count', type=int, default=1,
		help='the number of VMs to launch')

	args = parser.parse_args()

	if args.sisyphus:
		assert args.workflow, 'need a --workflow name to launch sisyphus workers'
		launch_sisyphus_workers(args.name_prefix, args.workflow, args.index, args.count)
	else:
		# TODO(jerry): + metadata?
		launch_fireworkers(args.name_prefix, args.index, args.count)


if __name__ == '__main__':
	main()
