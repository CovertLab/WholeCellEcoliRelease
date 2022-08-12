# Script to run simulations used to generate data for figures in the operon
# paper.

lpad reset
make clean compile

## Set A - no operons, glucose minimal media
# Used for comparisons in Figures 1, 2, 3, 4, 5
DESC="SET A 8 gens 128 seeds operons off with glucose minimal media" \
OPERONS="off" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=128 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set B - operons v1, glucose minimal media
# Used for comparisons in Figures 2
DESC="SET B 8 gens 128 seeds operons v1 with glucose minimal media" \
OPERONS="v1" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=128 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set C - operons v2, glucose minimal media
# Used for comparisons in Figures 2, 3
DESC="SET C 8 gens 128 seeds operons v2 with glucose minimal media" \
OPERONS="v2" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=128 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set D - operons v3, glucose minimal media
# Used for comparisons in Figure 3
DESC="SET D 8 gens 128 seeds operons v3 with glucose minimal media" \
OPERONS="v3" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=128 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set E - operons on (final version), glucose minimal media
# Used for comparisons in Figures 1, 4, 5
DESC="SET E 8 gens 128 seeds operons on with glucose minimal media" \
OPERONS="on" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=128 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py


## Launch the fireworks created with fw_queue.py
# Uncomment one method - rlaunch is interactive, qlaunch is distributed
# rlaunch rapidfire
# qlaunch -r rapidfire --nlaunches infinite --sleep 5 --maxjobs_queue 100
