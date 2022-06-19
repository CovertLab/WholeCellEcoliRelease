#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

import os

from six.moves import input

import wholecell.utils.filepath as fp


def default_input(prompt, default):
	'''
	Prompts the user for input with a default value that is returned if no
	input is received.

	Inputs:
		prompt (str) - prompt to display to user
		default (str) - default value accepted if no input received

	Returns str response from user or default if no input
	'''

	response = input("{} (default: {}): ".format(prompt, default))
	if response == "":
		print("  Using default: {}".format(default))
		response = default

	return response


def main():
	home = os.environ["HOME"]
	user = os.environ["USER"]

	print("\nEnter the following information for your FireWorks launchpad.  Hit return to accept defaults.\n")

	logdir_launchpad = default_input("Launchpad logging directory", os.path.join(home, "fw", "logs", "launchpad"))
	logdir_qadapter = default_input("Queue adapter logging directory", os.path.join(home, "fw", "logs", "qadapter"))

	db_host = default_input("\nDatabase host, e.g. ds257500.mlab.com", "localhost")
	db_port = default_input("Database port", "27017")

	db_name = default_input("\nDatabase name", user)
	db_username = default_input("Database username (empty if no login needed)", "")
	db_password = ""
	if db_username:
		db_password = default_input("Database password (stored in plaintext, unfortunately", "")

	wcecoli_path = fp.ROOT_PATH

	template_my_launchpad = os.path.join(wcecoli_path, "wholecell", "fireworks", "templates", "my_launchpad.yaml")
	my_launchpad = os.path.join(wcecoli_path, "my_launchpad.yaml")

	fp.makedirs(logdir_launchpad)
	fp.makedirs(logdir_qadapter)

	with open(template_my_launchpad, "r") as f:
		t = f.read()

	my_launchpad_text = t.format(
		LOGDIR_LAUNCHPAD=logdir_launchpad,
		DB_HOST=db_host,
		DB_NAME=db_name,
		DB_USERNAME=db_username or 'null',
		DB_PASSWORD=db_password or 'null',
		DB_PORT=db_port,
		)

	# TODO(jerry): Use a yaml-writer to get special cases right.
	with open(my_launchpad, "w") as f:
		f.write(my_launchpad_text)

	template_my_qadapter = os.path.join(wcecoli_path, "wholecell", "fireworks", "templates", "my_qadapter.yaml")
	my_qadapter = os.path.join(wcecoli_path, "my_qadapter.yaml")

	with open(template_my_qadapter, "r") as f:
		t = f.read()

	my_qadapter_text = t.format(
		LOGDIR_QADAPTER=logdir_qadapter,
		LAUNCHPAD_PATH=my_launchpad,
		WCECOLI_PATH=wcecoli_path,
		)

	with open(my_qadapter, "w") as f:
		f.write(my_qadapter_text)

	print("")
	print("Created {} with the information provided.".format(my_launchpad))
	print("Created {} with the information provided.".format(my_qadapter))


if __name__ == "__main__":
	main()
