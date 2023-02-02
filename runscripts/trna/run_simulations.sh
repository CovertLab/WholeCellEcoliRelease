# Script to run simulations used to generate data shown in the figures
# of the paper:
# todo

lpad reset
make clean compile

########################################################################
# Parameter Calculator -- uncomment one of the following 2 options:

# Uncomment to use optimized tRNA synthetase kinetic parameters
# python runscripts/manual/runParca.py sims

# Uncomment to run optimization (generates new tRNA synthetase kinetic
# parameters)
# python runscripts/manual/runParca.py --optimize-trna-charging sims

########################################################################
# Simulations

# Updated Model, Optimized Kinetics
python runscripts/manual/runSim.py -v trna_synthetase_kinetics 0 0 -g 10 -i 15 sims

# Prior Model, Unlimited tRNA Aminoacylation
python runscripts/manual/runSim.py -v trna_synthetase_kinetics 1 1 -g 10 -i 10 sims

# Updated Model, Highest Measured tRNA Synthetase kcats
python runscripts/manual/runSim.py -v trna_synthetase_kinetics 2 2 -g 2 -i 50 sims

# Updated Model, HisRS kcat = 142 1/s (highest measured value)
python runscripts/manual/runSim.py -v trna_synthetase_kinetics 3 3 -g 10 -i 15 sims

# Updated Model, HisRS kcat = 183 1/s
python runscripts/manual/runSim.py -v trna_synthetase_kinetics 4 4 -g 10 -i 10 sims

# Updated Model, HisRS kcat = 223 1/s
python runscripts/manual/runSim.py -v trna_synthetase_kinetics 5 5 -g 10 -i 10 sims

# Updated Model, ArgRS kcat = 26 1/s (highest measured value)
python runscripts/manual/runSim.py -v trna_synthetase_kinetics 6 6 -g 20 -i 10 sims
