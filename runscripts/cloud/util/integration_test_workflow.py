#!/usr/bin/env python

from runscripts.cloud.util.workflow_cli import WorkflowCLI


class TestWorkflow(WorkflowCLI):
	"""A test workflow for integration and regression tests of the cloud
	workflow software, which runs workflow tasks in Docker containers.

	The run will never reach "WORKFLOW COMPLETE" since some tasks intentionally
	fail and another task will never run for lack of inputs from failed tasks.
	"""

	DEFAULT_TIMEOUT = 30  # in seconds

	def build(self, args):
		"""Build the workflow."""
		# `print()` tests that a blank line in the console output doesn't cause
		# a confusing Logs Viewer message like this:
		# {"python_logger":"dockerfiretask.lines","message":""}
		lines_filename = self.internal('lines.txt')
		code = (
			"with open('" + lines_filename + "', 'w') as f:\n"
			"  for i in range(10):\n"
			"    f.write('This is line {}\\n'.format(i))\n"
			"    if i == 4:\n"
			"        print()\n"
			"    print('line {}'.format(i))")
		self.add_task(
			name='lines',
			outputs=[lines_filename],
			command=['python', '-u', '-c', code])

		self.add_task(
			name='count',
			inputs=[lines_filename],
			command=['wc', lines_filename])

		# Expected:  "wc: /tmp/lines.txt: No such file or directory"
		# because this task spec didn't request the input file.
		error_no_such_file_out = self.internal('expected_no_such_file.txt')
		self.add_task(
			name='expected_no_such_file',
			inputs=[],
			outputs=['>' + error_no_such_file_out],
			command=['wc', lines_filename])

		# Expected:  "IndexError: tuple index out of range"
		index_error_out = self.internal('expected_index_out_of_range.txt')
		self.add_task(
			name='expected_index_out_of_range',
			outputs=['>' + index_error_out],
			command=['python', '-u', '-c', "()[1]"])

		# This task depends on error task outputs as a regression test to check
		# that the workflow engine doesn't keep retrying any of them.
		# Expected:  This task never runs since its inputs never arrive.
		self.add_task(
			name='expected_to_never_run',
			inputs=[error_no_such_file_out, index_error_out],
			command=['cat', error_no_such_file_out, index_error_out])

		# This Task writes files into an output dir to test access from inside
		# the Container to an output dir created by the Fireworker outside the
		# Container and vice versa.
		# Expected:  The Fireworker can read and delete the Task's output files
		# without error.
		output_dir = self.internal('output/dir/')
		code = (
			"for i in range(4):\n"
			"  name = '" + output_dir + "{}.txt'.format(i)\n"
			"  print('Wrote ' + name)\n"
			"  with open(name, 'w') as f:\n"
			"    f.write('This is file {}\\n'.format(i))\n")
		self.add_task(
			name='output_dir',
			inputs=(),
			outputs=(output_dir,),
			command=['python', '-u', '-c', code])

		# Download and append to a file written by a previous Task to test that
		# the Task has write access to its input files, e.g. they're not owned
		# by root. Writing the file back to GCS with the same name would make
		# it ambiguous which Task is responsible for writing it, which doesn't
		# fit the functional WF data flow model. Skirt that by printing the
		# appended file to stdout and uploading stdout to demonstrate success.
		code = (
			"fn = '" + output_dir + "1.txt'\n"
			"with open(fn, 'a') as f:\n"
			"  f.write('Appended to file 1\\n')\n"
			"with open(fn, 'r') as f:\n"
			"  print(f.read())\n")
		self.add_task(
			name='overwrite',
			inputs=(output_dir,),
			outputs=['>' + self.internal('overwrite.txt')],
			command=['python', '-u', '-c', code])

		# test a timeout
		code = (
			"from time import sleep\n"
			"for i in range(1, 100):\n"
			"  sleep(1)\n"
			"  print('{:3} seconds'.format(i))")
		self.add_task(
			name='expected_timeout',
			timeout=8,
			command=['python', '-u', '-c', code])

		# Push the task (in Docker) to run out of memory (OOM) with 8GiB of data.
		# Spread it across multiple allocations so Python doesn't just punt.
		code = (
			"data = [bytearray(1024 * 1024 * 1024) for _ in range(8)]\n"
			"print(len(data), len(data[0]), data[0][-1])\n"
		)
		self.add_task(
			name='expected_out_of_memory',
			timeout=15,
			command=['python', '-u', '-c', code])


if __name__ == '__main__':
	TestWorkflow().cli()
