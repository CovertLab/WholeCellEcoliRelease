#!/bin/sh
# Measure Python library performance with several choices of # BLAS threads.
# Run this from the wcEcoli directory.

run() {
    echo OPENBLAS_NUM_THREADS=$1
    OPENBLAS_NUM_THREADS=$1 python -m wholecell.tests.utils.test_library_performance
}

echo "pyenv version: \c"
pyenv version

echo "SHERLOCK=$SHERLOCK"

echo "Running test_library_performance"

run ""
run 1
run 4
run 8
run 16
