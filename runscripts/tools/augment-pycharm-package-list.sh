#!/bin/sh
# Add to PyCharm's list of packages that don't match their import names to
# avoid inspection warnings like this:
#    Package containing module 'Bio' is not listed in project requirements
#
# Run this then restart PyCharm each time you install a PyCharm release.
# See https://youtrack.jetbrains.com/issue/PY-27985
#
# This handles PyCharm Pro (python.jar) and PyCharm CE (python-ce.jar)
# installed on macOS by JetBrains Toolbox.
#
# TODO: Support other PyCharm installers and other OSs. But this approach doesn't
# work for PyCharm installed by snap since it's on a read-only filesystem.

set -eu

FILE_TO_PATCH=tools/packages

patch_packages() {  # Patch the $FILE_TO_PATCH list inside the $1/$2 jar file.
  if [ -d "$1" ]; then
    cd "$1"

    if [ -e "$2" ]; then
      echo Checking PyCharm python*.jar file "$2"...

      unzip -o "$2" $FILE_TO_PATCH

      if grep --quiet vivarium_cell $FILE_TO_PATCH; then
        echo "$2" // $FILE_TO_PATCH was already patched
      else
        echo Patching "$2" // $FILE_TO_PATCH
        echo "
arrow stochastic-arrow
Bio biopython
borealis borealis-fireworks
vivarium vivarium-core
vivarium_cell vivarium-cell" >> $FILE_TO_PATCH
        zip python.jar $FILE_TO_PATCH
        echo "===Remember to restart PyCharm==="
      fi

      rm $FILE_TO_PATCH
    fi
  fi
}

MAC_JB_TOOLBOX_PYCHARM="$HOME/Library/Application Support/JetBrains/Toolbox/apps/PyCharm-P/ch-0"
if [ -d "$MAC_JB_TOOLBOX_PYCHARM" ]; then
  cd "$MAC_JB_TOOLBOX_PYCHARM"
  cd "$(ls -t | head -1)"
  patch_packages PyCharm.app/Contents/plugins/python/lib python.jar
  patch_packages PyCharm.app/Contents/plugins/python/lib python-ce.jar
fi
