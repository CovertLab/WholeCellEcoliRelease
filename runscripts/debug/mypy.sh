#!/bin/sh
# Run the mypy type checker on the wcEcoli sources (except prototypes/, etc.).
# The mypy.ini file sets some configuration options.
#
# ASSUMES: The current working directory is your wcEcoli git directory.

mypy
