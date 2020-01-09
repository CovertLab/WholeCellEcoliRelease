#!/usr/bin/env python
"""Launch Google Compute Engine VMs."""

### TODO: Deprecated. Switch to the replacement module in borealis-fireworks. ###


from __future__ import absolute_import, division, print_function

import argparse
import os
from pprint import pprint
import re
import sys

if os.name == 'posix' and sys.version_info[0] < 3:
	import subprocess32 as subprocess
else:
	import subprocess

from typing import Any, Dict, List, Optional
from wholecell.utils import filepath as fp


MAX_VMS = 100  # don't launch more than this many GCE VMs at a time

SERVICE_ACCOUNT = '441871726775-compute@developer.gserviceaccount.com'

SCOPES = ('https://www.googleapis.com/auth/devstorage.read_only,'
		  'https://www.googleapis.com/auth/logging.write,'
		  'https://www.googleapis.com/auth/monitoring.write,'
		  'https://www.googleapis.com/auth/servicecontrol,'
		  'https://www.googleapis.com/auth/service.management.readonly,'
		  'https://www.googleapis.com/auth/trace.append')

verbose = False


def gcloud_get_config(section_property):
	# type: (str) -> str
	"""Get a "section/property" configuration property from the gcloud command
	line tool. Raise an exception if the property is not configured.
	"""
	return fp.run_cmd(['gcloud', 'config', 'get-value', str(section_property)])


def gcp_project():
	# type: () -> str
	"""Get the current Google Cloud Platform (GCP) project from gcloud."""
	return gcloud_get_config('core/project')


def gce_zone():
	# type: () -> str
	"""Get the current Google Compute Engine (GCE) zone from gcloud."""
	return gcloud_get_config('compute/zone')


def make_VM_names(name_prefix, base=0, count=1):
	# type: (str, int, int) -> List[str]
	"""Make a list of valid GCE VM names from the prefix, base number, and count."""
	sanitized = re.sub(
		r'[^-a-z0-9]+', '-',
		name_prefix.lower().replace('workflow', ''))
	names = ['{}-{}'.format(sanitized, i) for i in range(base, base + count)]
	return names


def launch_GCE_VMs(name_prefix, base=0, count=1, command_options=None,
		dry_run=False, **metadata):
	# type: (str, int, int, Optional[Dict[str, Any]], bool, **Any) -> None
	"""Parallel-launch VM instances on GCE.

	This provides default command options for `gcloud compute instances create`.
	The caller should at least set the `image-family` option and override
	default options as needed.

	This converts command_options and metadata to strings, then passes tokens
	to `gcloud` avoiding quoting problems, but the metadata fields get joined
	with commas so those keys and values must not contain commas.
	"""
	def clean(token):
		# type: (Any) -> str
		"""Clean the token of "=" and "," chars so it won't mess up in a
		"key=val,key=val" metadata string.
		"""
		return re.sub(r'[=,]+', '', str(token))

	assert name_prefix, 'the name_prefix must not be empty'
	assert 0 <= count < MAX_VMS, 'VM count ({}) must be in the range [0 .. {}]'.format(
		count, MAX_VMS)
	instance_names = make_VM_names(name_prefix, base, count)

	dry = 'Dry run for: ' if dry_run else ''
	vms = ('VMs' if count <= 0
		   else 'VM: ' + instance_names[0] if count == 1
		   else 'VMs: {} .. {}'.format(instance_names[0], instance_names[-1]))
	print('{}Launching {} Google Compute Engine {}'.format(dry, count, vms))
	if count <= 0:
		return

	project = gcp_project()
	options = {
		'project': project,
		'image-project': project,
		'zone': gce_zone(),
		'machine-type': 'n1-standard-1',
		'subnet': 'default',
		'network-tier': 'PREMIUM',
		'maintenance-policy': 'MIGRATE',
		'boot-disk-size': '200GB',
		'boot-disk-type': 'pd-standard',
		'service-account': SERVICE_ACCOUNT,
		'scopes': SCOPES}
	options.update(command_options or {})

	metadata_string = ','.join(
		'{}={}'.format(clean(k), clean(v)) for k, v in metadata.items())
	if metadata_string:
		options['metadata'] = metadata_string

	options_list = ['--{}={}'.format(clean(k), v) for k, v in options.items()]
	cmd_tokens = ['gcloud', 'compute', 'instances', 'create'
		] + instance_names + options_list

	if dry_run or verbose:
		pprint(cmd_tokens)

	if not dry_run:
		subprocess.call(cmd_tokens)  #, env=os.environ


def launch_sisyphus_workers(name_prefix, workflow, base=0, count=1,
		command_options=None, dry_run=False, **metadata):
	# type: (str, str, int, int, Optional[Dict[str, Any]], bool, **Any) -> None
	"""Parallel-launch Sisyphus worker VMs on GCE for the given workflow name."""
	assert workflow, 'a workflow name is required'

	options = {
		'image-family': 'sisyphus-worker',
		'description': 'sisyphus worker'}
	options.update(command_options or {})

	launch_GCE_VMs(name_prefix, base, count, options, dry_run=dry_run,
		workflow=workflow, **metadata)


def launch_fireworkers(name_prefix, base=0, count=1, command_options=None,
		dry_run=False, **metadata):
	# type: (str, int, int, Optional[Dict[str, Any]], bool, **Any) -> None
	"""Parallel-launch FireWorks worker VMs on GCE."""
	options = {
		'image-family': 'fireworker',
		'description': 'FireWorks worker'}
	options.update(command_options or {})

	launch_GCE_VMs(name_prefix, base, count, options, dry_run=dry_run, **metadata)


def main():
	parser = argparse.ArgumentParser(
		description='Launch Google Compute Engine VMs all at once.'
					' Python code can call the launch_...() functions instead'
					' of this command line interface.')
	parser.add_argument('name_prefix', metavar='NAME-PREFIX',
		help='The GCE VM name prefix for constructing VM names in the pattern'
			 ' {PREFIX}-{NUMBER}.')
	parser.add_argument('-b', '--base', type=int, default=0,
		help='The base number for the numbered VM names (default 0). Use this'
			 ' to launch additional VMs with unique names.')
	parser.add_argument('-c', '--count', type=int, default=1,
		help='The number of VMs to launch (default 1).')
	parser.add_argument('-d', '--dry-run', action='store_true', dest='dry_run',
		help='Dry run: Print the `gcloud` command, then exit.')
	parser.add_argument('-s', '--sisyphus', action='store_true',
		help='Launch Sisyphus workers (rather than Fireworkers).')
	parser.add_argument('-w', '--workflow',
		help='The workflow name to use with --sisyphus.')
	parser.add_argument('-m', '--metadata', metavar='KEY=VALUE', nargs='*',
		default=[],
		help='Custom metadata settings, e.g. "db=crick" to identify a FireWorks'
			 ' LaunchPad MongoDB database to the workers.')

	args = parser.parse_args()
	unpacked = [e.split('=', 2) + [''] for e in args.metadata]
	metadata = {e[0]: e[1] for e in unpacked}

	if args.sisyphus:
		assert args.workflow, 'need a --workflow name to launch sisyphus workers'
		launch_sisyphus_workers(args.name_prefix, args.workflow, args.base,
			args.count, dry_run=args.dry_run, **metadata)
	else:
		launch_fireworkers(args.name_prefix, args.base, args.count,
			dry_run=args.dry_run, **metadata)


if __name__ == '__main__':
	main()
