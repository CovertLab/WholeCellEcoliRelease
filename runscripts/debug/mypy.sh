#!/bin/sh
# Run the mypy type checker on the wcEcoli sources (except prototypes/, etc.)
# The mypy.ini file sets some configuration options.
# So far it only runs a check for Python 2.7 check -- Python 3 checks to come.
# A key goal is to find problems for Python 3 in unicode vs. bytes and
# modernize code changes such as range() becoming an iterator.
#
# mypy now runs clean (thanks to skipping some cases) so we should run it
# before making a Pull Request and in CI.
#
# It's also smart to run PyCharm inspections on the code.
#
# ASSUMES: The current working directory is your wcEcoli git directory.
#
# ASSUMES: You have a `mypy` virtualenv, created like this:
#   pyenv install 3.8.2
#   pyenv shell 3.8.2
#   pyenv virtualenv mypy
#   pyenv shell mypy
#   pip install mypy
#   pyenv shell --unset

# Run in the Python 3 virtualenv `mypy` with access to the current virtualenv,
# which is typically the Python 2.7 virtualenv `wcEcoli2`.
export PYENV_VERSION="mypy:$(pyenv version-name)"

mypy --py2
