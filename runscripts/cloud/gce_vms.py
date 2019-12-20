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


def launch_GCE_VMs(instance_names, command_options=None, **metadata):
	# type: (Iterable[str], Optional[Dict[str, Any]], **Any) -> None
	"""Parallel-launch VMs on GCE. This function provides default options for
	project, zone, etc. Do set at least image-family in command_options and
	override defaults as needed."""
	assert not isinstance(instance_names, str), (
		'instance_names must be a collection of names')
	instance_names = list(instance_names)[:MAX_VMS]
	for name in instance_names:
		assert re.match(r'[-a-z0-9]+$', name), '"{}" is not a valid GCE VM name'.format(name)

	print('\nLaunching {} Google Compute Engine VM(s).'.format(len(instance_names)))
	if len(instance_names) <= 0:
		return

	if not command_options:
		command_options = {}

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
	options.update(command_options)

	metadata_string = ','.join('{}={}'.format(k, v) for k, v in metadata.items())
	if metadata_string:
		options['metadata'] = metadata_string

	options_list = ['--{}={}'.format(k, v) for k, v in options.items()]
	cmd_tokens = ['gcloud', 'compute', 'instances', 'create'
		] + instance_names + options_list
	subprocess.call(cmd_tokens, env=os.environ)


def make_VM_names(prefix, base=0, count=1):
	# type: (str, int, int) -> List[str]
	"""Make a list of valid GCE VM names from the prefix, base number, and count."""
	sanitized = re.sub(r'[^-a-z0-9]+', '-', prefix.lower().replace('workflow', ''))
	names = ['{}-{}'.format(sanitized, i) for i in range(base, base + count)]
	return names


def launch_sisyphus_workers(instance_names, workflow):
	# type: (Iterable[str], str) -> None
	"""Parallel-launch Sisyphus worker VMs on GCE for the given workflow name."""
	assert workflow, "a workflow name is required"

	options = {'image-family': 'sisyphus-worker'}
	launch_GCE_VMs(instance_names, options, workflow=workflow, description='sisyphus worker')


def launch_fireworkers(instance_names, **metadata):
	# type: (Iterable[str], **Any) -> None
	"""Parallel-launch FireWorks worker VMs on GCE."""
	options = {'image-family': 'fireworker'}
	launch_GCE_VMs(instance_names, options, description='fireworker', **metadata)


def main():
	parser = argparse.ArgumentParser(
		description='Launch Google Compute Engine VMs.')
	parser.add_argument('--sisyphus', action='store_true',
		help='Launch Sisyphus workers (rather than Fireworkers)')
	parser.add_argument('-w', '--workflow',
		help='the workflow name; required for --sisyphus')
	parser.add_argument('names', nargs='+', help='the GCE VM names to launch')

	args = parser.parse_args()

	if args.sisyphus:
		assert args.workflow, 'need a --workflow name to launch sisyphus workers'
		launch_sisyphus_workers(args.names, args.workflow)
	else:
		launch_fireworkers(args.names)  # TODO(jerry): optional metadata?


if __name__ == '__main__':
	main()

