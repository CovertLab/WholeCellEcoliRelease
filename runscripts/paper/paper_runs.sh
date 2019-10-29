# Script to run simulations used to generate data for figures in the paper
# Data will be uncompressed to allow running of new analysis scripts
# Requires setup outlined in requirements.txt and in README.md
# Suggestion: run one set at a time by commenting out all other fw_queue.py lines

lpad reset
make clean compile

## Set A - basal condition
# Used for figure 4
DESC="SET A 32 gens 8 seeds basal with growth noise and D period" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=32 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
python runscripts/fw_queue.py

## Set B - nutrient shift from minimal to minimal + AA
# Used for figure 1 and S1
DESC="SET B 9 gens 8 seeds shift to plus AA without growth noise with D period" \
VARIANT="nutrientTimeSeries" FIRST_VARIANT_INDEX=2 LAST_VARIANT_INDEX=2 \
SINGLE_DAUGHTERS=1 N_GENS=9 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=0 D_PERIOD_DIVISION=1 \
python runscripts/fw_queue.py

## Set C - 3 growth rates from different conditions
# Used for figure 2 and S2
DESC="SET C 4 gens 256 seeds 3 conditions with growth noise and D period" \
VARIANT="condition" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=2 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=256 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
RUN_AGGREGATE_ANALYSIS=0 WC_LENGTHSEC=14400 \
python runscripts/fw_queue.py

## Set D - changes to RNAP and ribosome expression
# Used for figure 2

# This subset of simulations is redundant with the first variant of Set C
DESC="SET D1 4 gens 256 seeds" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=256 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
DISABLE_RIBOSOME_CAPACITY_FITTING=0 DISABLE_RNAPOLY_CAPACITY_FITTING=0 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

DESC="SET D2 4 gens 256 seeds, unfit ribosome expression" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=256 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=0 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

DESC="SET D3 4 gens 256 seeds, unfit rna poly expression" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=256 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
DISABLE_RIBOSOME_CAPACITY_FITTING=0 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

DESC="SET D4 4 gens 256 seeds, unfit ribosome and rna poly expression, adjusted RNase expression" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=256 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
ADJUST_RNASE_EXPRESSION=1 RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

# This subset of simulations is redundant with the first variant of Set L
DESC="SET D5 4 gens 256 seeds, unfit ribosome and rna poly expression" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=256 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

## Set E - metabolism objective weighting
# Used for figure 3 and S3
DESC="SET E 8 gens 8 seeds 10 metabolism weighting values" \
VARIANT="metabolism_kinetic_objective_weight" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=9 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=0 D_PERIOD_DIVISION=0 \
python runscripts/fw_queue.py

## Set F - pabB expression
# Used in figure S4
DESC="SET F 16 gens 8 seeds 9 pabB expression values" \
VARIANT=subgen_expression FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=8 \
SINGLE_DAUGHTERS=1 N_GENS=16 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

## Set G - no RNA and protein parameter adjustments
DESC="SET G 32 gens 8 seeds basal without parameter adjustments" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=32 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
ADJUST_RNA_AND_PROTEIN_PARAMETERS=0 DISABLE_MEASURED_PROTEIN_DEG=1 \
python runscripts/fw_queue.py

## Set H - setup for metabolism_kinetic_objective_interactions variant
# Use python models/ecoli/analysis/variant/flux_sensitivity.py to identify
# reactions for factorial design of experiments
DESC="SET H 1 gen flux sensitivity" \
VARIANT="flux_sensitivity" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=1 N_INIT_SIMS=1 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=0 D_PERIOD_DIVISION=0 \
RUN_AGGREGATE_ANALYSIS=0 WC_LENGTHSEC=10 \
python runscripts/fw_queue.py

## Set I - factorial design experiments for disabled kinetic constraints
DESC="SET I kinetic constraint factorial design" \
VARIANT="kinetic_constraints_factorial_experiments" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=256 \
SINGLE_DAUGHTERS=1 N_GENS=1 N_INIT_SIMS=4 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=0 D_PERIOD_DIVISION=0 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

## Set J - changes to input data sets for parameters described below:
# Used in figures 2 and S2
# RNA mass per cell
DESC="SET J01 alternate RNA mass per cell" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
ALTERNATE_MASS_FRACTION_RNA=1 RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

# mRNA mass fraction
DESC="SET J02 alternate mRNA mass per cell" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
ALTERNATE_MASS_FRACTION_MRNA=1 RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

# RNA-sequencing (This set is redundant with SET D1, but fewer seeds)
DESC="SET J03 alternate RNA seq" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

# Translation efficiency
DESC="SET J04 alternate translation efficiency" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
ALTERNATE_TRANSLATION_EFFICIENCY=1 RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

# R-protein half lives
DESC="SET J05 alternate R-protein half lives" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
ALTERNATE_R_PROTEIN_DEGRADATION=1 RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

# Protein mass per cell
DESC="SET J06 alternate protein mass per cell" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
ALTERNATE_MASS_FRACTION_PROTEIN=1 RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

# Transcription elongation rate
DESC="SET J07 variable transcription elongation rate" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
VARIABLE_ELONGATION_TRANSCRIPTION=1 RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

# RNA polymerase active fraction
DESC="SET J08 max RNA polymerase activity" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
MAX_RNAP_ACTIVITY=1 RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

# mRNA half life
DESC="SET J09 mRNA half life fitting" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
MRNA_HALF_LIFE_FITTING=1 RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

# Translation elongation rate
DESC="SET J10 variable translation elongation rate" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
VARIABLE_ELONGATION_TRANSLATION=1 RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

# Ribosome active fraction
DESC="SET J11 alternate ribosome activity" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
ALTERNATE_RIBOSOME_ACTIVITY=1 RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

## Set K - parameter sensitivity analysis
runscripts/paper/sensitivity.sh "SET_K_sensitivity" 0 20000

## Set L - comparison of RNA/protein ratios with Scott et al. before changes
# to RNAP and ribosome expression
# Used for Figure S2
DESC="SET L 4 gens 256 seeds 3 conditions unfit ribosome and rna poly expression" \
VARIANT="condition" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=2 \
SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=256 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \
DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

## Set M - larger sims to compare kinetic constraints
# Used for Figure 3
DESC="SET M1 old kinetic constraints" \
VARIANT="kinetic_constraints_factorial_experiments" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=0 D_PERIOD_DIVISION=0 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

DESC="SET M2 new kinetic constraints" \
VARIANT="kinetic_constraints_factorial_experiments" FIRST_VARIANT_INDEX=47 LAST_VARIANT_INDEX=47 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=8 \
MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=0 D_PERIOD_DIVISION=0 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fw_queue.py

## Launch the fireworks created with fw_queue.py
# Uncomment one method - rlaunch is interactive, qlaunch is distributed
# rlaunch rapidfire
# qlaunch -r rapidfire --nlaunches infinite --sleep 5 --maxjobs_queue 100
