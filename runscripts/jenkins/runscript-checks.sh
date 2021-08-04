#! /usr/bin/env bash

# Test some runscripts to make sure they do not fall out of step with other updates
# Remove any of these if they are no longer useful and causing build failures

set -e

# These should come before parameter_search.py to use the correct default directory
runscripts/debug/metabolism.py --validation -1
runscripts/debug/charging.py --validation -1

runscripts/manual/parameter_search.py --method quick_example
