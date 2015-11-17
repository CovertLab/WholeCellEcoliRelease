#!/usr/bin/env python

import os
from jinja2 import Template

def main():

	home = os.environ["HOME"]

	print "The following values are for your primary launchpad, probably the one you made first on mongolab.com."
	logdir_launchpad = raw_input("Launchpad logging directory (e.g. %s): " % os.path.join(home, "fw", "logs", "launchpad"))
	db_host = raw_input("Database host (e.g., x.mongolab.com): ")
	db_name = raw_input("Database name (e.g., wc_ecoli): ")
	db_username = raw_input("Database username (e.g., fireworks): ")
	db_password = raw_input("Database password (stored in plaintext, unfortunately): ")
	db_port = raw_input("Database port: ")

	logdir_qadapter = raw_input("Queue adapter logging directory (e.g. %s): " % os.path.join(home, "fw", "logs", "qadapter"))
	wcecoli_path = raw_input("wcEcoli path (e.g., %s): " % os.path.join(home, "wcEcoli"))

	template_my_launchpad = os.path.join(wcecoli_path, "wholecell", "fireworks", "templates", "my_launchpad.yaml")
	my_launchpad = os.path.join(wcecoli_path, "my_launchpad.yaml")

	if not os.path.exists(logdir_launchpad):
		os.makedirs(logdir_launchpad)
	if not os.path.exists(logdir_qadapter):
		os.makedirs(logdir_qadapter)

	h = open(template_my_launchpad, "r")
	t = Template(h.read())
	h.close()

	my_launchpad_text = t.render({
		"LOGDIR_LAUNCHPAD": logdir_launchpad,
		"DB_HOST": db_host,
		"DB_NAME": db_name,
		"DB_USERNAME": db_username,
		"DB_PASSWORD": db_password,
		"DB_PORT": db_port,
		})

	h = open(my_launchpad, "w")
	h.write(my_launchpad_text)
	h.close()

	template_my_qadapter = os.path.join(wcecoli_path, "wholecell", "fireworks", "templates", "my_qadapter.yaml")
	my_qadapter = os.path.join(wcecoli_path, "my_qadapter.yaml")

	h = open(template_my_qadapter, "r")
	t = Template(h.read())
	h.close()

	my_qadapter_text = t.render({
		"LOGDIR_QADAPTER": logdir_qadapter,
		"LAUNCHPAD_PATH": my_launchpad,
		})

	h = open(my_qadapter, "w")
	h.write(my_qadapter_text)
	h.close()
	
	print ""
	print "Created %s with the information provided." % my_launchpad
	print "Created %s with the information provided." % my_qadapter

	print ""
	print "The following values are for your secondary launchpad, probably the one you made second on mongolab.com"
	db_host = raw_input("Secondary launchpad database host (e.g., x.mongolab.com): ")
	db_name = raw_input("Secondary launchpad database name (e.g., wc_ecoli_2): ")
	db_username = raw_input("Secondary launchpad database username (e.g., fireworks): ")
	db_password = raw_input("Secondary launchpad database password (stored in plaintext, unfortunately): ")
	db_port = raw_input("Secondary launchpad database port: ")

	template_my_launchpad = os.path.join(wcecoli_path, "wholecell", "fireworks", "templates", "my_launchpad.yaml")
	my_launchpad = os.path.join(wcecoli_path, "my_launchpad_2.yaml")

	if not os.path.exists(logdir_launchpad):
		os.makedirs(logdir_launchpad)
	if not os.path.exists(logdir_qadapter):
		os.makedirs(logdir_qadapter)

	h = open(template_my_launchpad, "r")
	t = Template(h.read())
	h.close()

	my_launchpad_text = t.render({
		"LOGDIR_LAUNCHPAD": logdir_launchpad,
		"DB_HOST": db_host,
		"DB_NAME": db_name,
		"DB_USERNAME": db_username,
		"DB_PASSWORD": db_password,
		"DB_PORT": db_port,
		})

	h = open(my_launchpad, "w")
	h.write(my_launchpad_text)
	h.close()

	template_my_qadapter = os.path.join(wcecoli_path, "wholecell", "fireworks", "templates", "my_qadapter.yaml")
	my_qadapter = os.path.join(wcecoli_path, "my_qadapter_2.yaml")

	h = open(template_my_qadapter, "r")
	t = Template(h.read())
	h.close()

	my_qadapter_text = t.render({
		"LOGDIR_QADAPTER": logdir_qadapter,
		"LAUNCHPAD_PATH": my_launchpad,
		})

	h = open(my_qadapter, "w")
	h.write(my_qadapter_text)
	h.close()
	
	print ""
	print "Created %s with the information provided." % my_launchpad
	print "Created %s with the information provided." % my_qadapter

if __name__ == "__main__":
	main()
