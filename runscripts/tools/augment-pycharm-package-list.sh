#!/bin/sh
# Add to PyCharm's list of packages that don't match their import names to
# avoid inspection warnings like this:
#    Package containing module 'Bio' is not listed in project requirements
#
# Run this then restart PyCharm each time you install a PyCharm release.
# See https://youtrack.jetbrains.com/issue/PY-27985
#
# This handles PyCharm Pro installed on macOS by JetBrains Toolbox.
# TODO: Support PyCharm CE (python-ce.jar), other installers, and other OSs.

set -eu

cd ~/Library/Application\ Support/JetBrains/Toolbox/apps/PyCharm-P/ch-0
cd $(ls -t | head -1)
cd PyCharm.app/Contents/plugins/python/lib

unzip python.jar tools/packages
echo "arrow stochastic-arrow
Bio biopython
borealis borealis-fireworks
vivarium wholecell-vivarium" >> tools/packages
zip python.jar tools/packages
