#!/bin/sh
# Run the mypy type checker on the wcEcoli sources (except prototypes/, etc.).
# This type-checks with the current pyenv Python -- 2.7 or 3.x.
# The mypy.ini file sets some configuration options.
# (It's also useful to run PyCharm inspections on the code.)
#
# ASSUMES: The current working directory is your wcEcoli git directory.
#
# ASSUMES: You have a `mypy` virtualenv created using Python 3.5+:
#   pyenv install 3.8.3
#   pyenv shell 3.8.3
#   pyenv virtualenv mypy
#   pyenv shell mypy
#   pip install mypy==0.770
#   pyenv shell --unset

PYTHON_EXECUTABLE="$(pyenv which python)"
PYTHON_VERSION=$(python -c "import sys; print('{}.{}'.format(*sys.version_info[:2]))")
PYENV_NAME=$(pyenv version-name)

# Run mypy within the Python 3 virtualenv `mypy` with access to the current virtualenv.
export PYENV_VERSION="mypy:$PYENV_NAME"
echo --- Running mypy for Python "$PYTHON_VERSION" in "$PYENV_NAME" ---
mypy --python-version "$PYTHON_VERSION" --python-executable "$PYTHON_EXECUTABLE"
