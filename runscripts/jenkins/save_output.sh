#! /bin/bash

# Script to move completed Jenkins runs to an output directory.
# Also adds the output to a tar archive in the output directory 3x per month
# to prevent purging of the archive and loss of data.

# Usage: ./save_output.sh src_dir dest_dir
#   src_dir: will save all directories inside and add to tar if date has not been saved
#   dest_dir: directory to move sims to and add to the existing tar archive inside

set -e

src_dir=$1
dest_dir=$2

# Check correct number of inputs
if [ $# -ne 2 ]; then
    echo "Must supply src and dest dir"
    exit
fi

# Make dest_dir an absolute path for when directory is changed
dest_dir=$(realpath $dest_dir)

# Ensure dest_dir exists
mkdir -p $dest_dir

# Get date string to only save in tar every 11 days (3x per month)
# (eg 20041 for all dates between 4/11/20 - 4/21/20)
date_str="$(date +%y%m)$((10#$(date +%d) / 11))"

# Files for tar backup
label=$(basename $dest_dir)
date_file="${dest_dir}/${label}_tar_dates.txt"
tar_file="${dest_dir}/${label}_sims.tar"

# Ensure date_file exists if a new directory
touch $date_file

# If current date string does not have a sim saved in tar,
# then add the output from this run to the archive
cd $src_dir
if [ -z "$(grep $date_str $date_file)" ]; then
    # If running with slurm, only add to the tar file if more than an hour
    # of execution time is left to prevent possible interruption during tar
    if [ -z "$SLURM_JOBID" ] || [ $(squeue -j $SLURM_JOBID -h -o '%L' | wc -c) -gt 6 ]; then
        echo "Adding sims to $tar_file"
        echo $date_str >> $date_file
        tar rf $tar_file *
    fi
fi

mv * $dest_dir
