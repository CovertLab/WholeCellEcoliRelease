#! /bin/bash

# Script to create an archive of completed sims similar to save_output.sh.
# Need to run from a directory containing sim output directories to add to
# the archive. This will add, at most, 3 sims per month to the archive.

set -e

# Files for tar backup
label=$(basename $(pwd))
date_file="${label}_tar_dates.txt"
tar_file="${label}_sims.tar"

# Ensure date file exists if a new directory
touch $date_file

# Find all directories in current directory starting with timestamp YYYYMMDD
# to add to tar archive
for src_dir in `find * -maxdepth 0 -type d -regextype sed -regex '^[0-9]\{8\}.*'`; do
    # Get date string to only save in tar every 11 days (3x per month)
    # (eg 20041 for all dates between 4/11/20 - 4/21/20)
    src_date=${src_dir::8}
    date_str="${src_date:2:4}$((10#${src_date:6:2} / 11))"

    # If current month does not have a sim saved in tar,
    # then add the output from all directories matching the current date
    # for cases like optional features that have multiple sims on a date.
    if [ -z "$(grep $date_str $date_file)" ]; then
        echo "Adding sims from $src_date to $tar_file"
        echo $date_str >> $date_file
        tar rf $tar_file ${src_date}*
    fi
done
