#!/bin/sh
# This script gets information about fireworks launched through qlaunch
# It will display the out or error file for the given fw_id
# Run from wcEcoli directory where block directories are stored
#
# Arguments:
#   - first position: "out" or "error" - select based on desired file
#   - second position: fw_id - ID of firework to get info about
#
# Example use:
#   ./runscripts/fireworks/fw_info.sh out 1 - gets output info for firework with ID of 1
#
# Notes:
#   - For multiple files with same ID (eg reran firework or reset launchpad), this will display each one
#   - TIP: add aliases out="./runscripts/fireworks/fw_info.sh out" and err="./runscripts/fireworks/fw_info.sh error"
#     to your .bash_profile so you can call with just `out 1` or `err 1`

fw="\"fw_id\": $2,"
for result in `grep -rl "$fw" block*`; do
    parent=$(dirname "$result")
    for f in `find "$parent" -name "*.$1"`; do
	echo $f
	less $f
    done
done

if [[ $result == '' ]]; then
    echo "Could not find fw_id $2"
fi
