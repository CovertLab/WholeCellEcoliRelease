#!/bin/bash

if [ -z "$OMPI_COMM_WORLD_RANK" ] || [ "$OMPI_COMM_WORLD_RANK" -eq 0 ]; then
	nosetests "$@" --with-xunit --xunit-file="nosetests.0.xml"
else
	nosetests "$@" --with-xunit --xunit-file="nosetests.$OMPI_COMM_WORLD_RANK.xml" #> /dev/null 2>&1
fi
