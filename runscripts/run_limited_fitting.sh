
# echo Submitting job 1 - fit all

# DISABLE_RIBOSOME_CAPACITY_FITTING=0 DISABLE_RNAPOLY_CAPACITY_FITTING=0 \
# PARALLEL_FITTER=1 \
# SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=1 DESC="Both fit" \
# python runscripts/fw_queue.py

# echo Submitting job 2 - fit only ribosomes

# DISABLE_RIBOSOME_CAPACITY_FITTING=0 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
# PARALLEL_FITTER=1 \
# SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=1 DESC="Ribosomes fit" \
# python runscripts/fw_queue.py

# echo Submitting job 3 - fit only RNA polymerases

# DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=0 \
# PARALLEL_FITTER=1 \
# SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=1 DESC="RNApoly fit" \
# python runscripts/fw_queue.py

echo Submitting job 4 - fit neither ribosomes nor RNA polymerases

DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
PARALLEL_FITTER=1 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=1 DESC="Neither fit" \
python runscripts/fw_queue.py
