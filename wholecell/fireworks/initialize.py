#!/usr/bin/env python

import os

def default_input(prompt, default):
	'''
	Prompts the user for input with a default value that is returned if no
	input is received.

	Inputs:
		prompt (str) - prompt to display to user
		default (str) - default value accepted if no input received

	Returns str response from user or default if no input
	'''

	response = raw_input("{} (e.g., {}): ".format(prompt, default))
	if response == "":
		print "Using default: {}".format(default)
		response = default

	return response

def main():

	home = os.environ["HOME"]

	print "Enter the following information for your launchpad.  Hitting return on directories will accept the default value."
	logdir_launchpad = default_input("Launchpad logging directory", os.path.join(home, "fw", "logs", "launchpad"))
	db_host = raw_input("Database host (e.g., x.mlab.com): ")
	db_name = raw_input("Database name (e.g., wc_ecoli): ")
	db_username = raw_input("Database username (e.g., fireworks): ")
	db_password = raw_input("Database password (stored in plaintext, unfortunately): ")
	db_port = raw_input("Database port: ")

	logdir_qadapter = default_input("Queue adapter logging directory", os.path.join(home, "fw", "logs", "qadapter"))
	wcecoli_path = default_input("wcEcoli path", os.path.join(home, "wcEcoli"))

	template_my_launchpad = os.path.join(wcecoli_path, "wholecell", "fireworks", "templates", "my_launchpad.yaml")
	my_launchpad = os.path.join(wcecoli_path, "my_launchpad.yaml")

	if not os.path.exists(logdir_launchpad):
		os.makedirs(logdir_launchpad)
	if not os.path.exists(logdir_qadapter):
		os.makedirs(logdir_qadapter)

	with open(template_my_launchpad, "r") as f:
		t = f.read()

	my_launchpad_text = t.format(
		LOGDIR_LAUNCHPAD=logdir_launchpad,
		DB_HOST=db_host,
		DB_NAME=db_name,
		DB_USERNAME=db_username,
		DB_PASSWORD=db_password,
		DB_PORT=db_port,
		)

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

	print ""
	print "Created {} with the information provided.".format(my_launchpad)
	print "Created {} with the information provided.".format(my_qadapter)

if __name__ == "__main__":
	main()
