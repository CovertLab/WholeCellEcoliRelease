#! /bin/bash
# Script to run many variants in parallel for parameter sensitivity analysis.
# Fireworks is too slow and checks too many task dependencies for this scale.
# Also removes variant sim_data objects to save space since affected parameters
# can be recalculated in the analysis script.

# Note: if parallel is not installed, simulations will be run sequentially,
#   to install on Linux: sudo apt install parallel

# Usage (from wcEcoli home directory):
#   runscripts/paper/sensitivity.sh [output dir] [start variant] [final variant]

set -u

# Required arguments
out_dir=$1
start_var=$2
end_var=$3

# Run simulation and remove sim_data to save space
# Args: output directory, log directory, variant to run
function simulation {
	variant="param_sensitivity"
	python runscripts/manual/runSim.py $1 --length_sec 10 --variant $variant $3 $3 > $2/$3.log 2>&1
	rm out/$1/${variant}_$(printf "%06d" $3)/kb/simData_Modified.cPickle
}
export -f simulation

# Create sim_data
python runscripts/manual/runParca.py $out_dir

# Create log directory
log_dir=out/$out_dir/log
if [ ! -e $log_dir ]; then
	mkdir $log_dir
fi

# Run simulation variants
## Use parallel if it exists
if [ -n "$(type -t parallel)" ]; then
    echo "$(date): Running simulations in parallel"
    seq $start_var $end_var | parallel simulation $out_dir $log_dir
## Otherwise run sequentially
else
    echo "$(date): Running simulations sequentially...try installing 'parallel' to speed up simulations"
    for i in $(seq $start_var $end_var); do
		simulation $out_dir $log_dir $i
    done
fi

# Run analysis script
echo "$(date): Running analysis"
python models/ecoli/analysis/variant/param_sensitivity.py $out_dir > $log_dir/analysis.log
