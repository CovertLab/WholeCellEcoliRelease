#!/bin/bash

INPUT="$1"
OUTPUT="$2"

echo $INPUT
echo $OUTPUT
cat $INPUT | sed 's/"""/"/g' | sed 's/""/"/g' | sed 's/"\[/\[/g' | sed 's/\]"/\]/g' | sed 's/"{/{/g' | sed 's/}"/}/g' | sed 's/1\/s/1 \/ units\.s/g' | sed 's/ (uM)//g' > $OUTPUT
