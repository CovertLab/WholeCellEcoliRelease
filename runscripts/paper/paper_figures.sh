# Script to run analysis to generate published figures
# When reproducing, /scratch/PI/mcovert/wc_ecoli/paper will be the out directory
# Simulations in runscripts/paper/paper_runs.sh must be run to generate the data before performing analysis


## Figure 1
#Bottom Simulations Left to Right, Top to bottom:
 python models/ecoli/analysis/multigen/growth_dynamics_grid.py \
 /scratch/PI/mcovert/wc_ecoli/paper/SET_B/20180822.181305.660233__SET_B_9_gens_8_seeds_shift_to_plus_AA_without_growth_noise_with_D_period/ \
 --variant_index 2 --seed 7


## Figure 2
# Panel a
python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_D/20180815.085246.549352__SET_D1_4_gens_256_seeds/

python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_D/20180815.181420.368837__SET_D2_4_gens_256_seeds,_unfit_ribosome_expression/

python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_D/20180815.235510.655007__SET_D3_4_gens_256_seeds,_unfit_rna_poly_expression/

python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_D/20180817.001302.336012__SET_D4_4_gens_256_seeds,_unfit_ribosome_and_rna_poly_expression/

# Panel b
python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20180810.094634.563490__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/ \
--variant_index 1

python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20180810.094634.563490__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/ \
--variant_index 2

# Panel c
python models/ecoli/analysis/variant/growth_condition_comparison_validation.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20180810.094634.563490__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/

# Panel d
python models/ecoli/analysis/variant/adder_sizer.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20180810.094634.563490__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/


## Figure 3
# panel A and D
python models/ecoli/analysis/variant/metabolism_kinetic_objective_weight.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20180813.102120.814493__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/

# panel B
python models/ecoli/analysis/single/massFractionSummary.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20180813.102120.814493__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 3 --generation 0 --seed 0

# panel C
python models/ecoli/analysis/single/massFractionSummary.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20180813.102120.814493__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 6 --generation 0 --seed 0

# panel E
python models/ecoli/analysis/cohort/kinetics_flux_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20180813.102120.814493__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 0

# panel F
python models/ecoli/analysis/cohort/kinetics_flux_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20180813.102120.814493__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 3

# panel G
python models/ecoli/analysis/cohort/centralCarbonMetabolismScatter.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20180813.102120.814493__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 3


## Figure 4
# panel A
python models/ecoli/analysis/multigen/proteinCountsValidation.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_A/20180809.155239.564816__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/ \
--seed 0

# panel B
python models/ecoli/analysis/cohort/expression_dynamics.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_A/20180809.155239.564816__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/

# panel C
python models/ecoli/analysis/multigen/subgenerationalTranscription.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_A/20180809.155239.564816__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/ \
--seed 0

# panel D
python models/ecoli/analysis/multigen/subgenerationalTranscription.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_A/20180809.155239.564816__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/ \
--seed 0

# panel E
python models/ecoli/analysis/multigen/functionalUnits_new.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_A/20180809.155239.564816__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/ \
--seed 0

# panel F
python models/ecoli/analysis/multigen/mene_limitations.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_A/20180809.155239.564816__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/ \
--seed 0


## Figure 5
# panel A: made by JC in repo EcoliFoldChanges

# panel B: made by JC from his experimental results


## Supplemental figure 1
python models/ecoli/analysis/multigen/growth_dynamics_panel.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_B/20180822.181305.660233__SET_B_9_gens_8_seeds_shift_to_plus_AA_without_growth_noise_with_D_period/ \
--variant_index 2 --seed 7

python models/ecoli/analysis/multigen/environmental_shift_fluxes.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_B/20180822.181305.660233__SET_B_9_gens_8_seeds_shift_to_plus_AA_without_growth_noise_with_D_period/ \
--seed 7


## Supplemental figure 2
# Made by script for panel d, Figure 2


## Supplemental figure 3
# panel A
python models/ecoli/analysis/single/massFractionSummary.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20180813.102120.814493__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 0 --generation 0 --seed 0

python models/ecoli/analysis/cohort/kinetics_flux_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20180813.102120.814493__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 0

python models/ecoli/analysis/cohort/centralCarbonMetabolismScatter.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20180813.102120.814493__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 0

# panel B
python models/ecoli/analysis/single/massFractionSummary.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20180813.102120.814493__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 3 --generation 0 --seed 0

python models/ecoli/analysis/cohort/kinetics_flux_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20180813.102120.814493__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 3

python models/ecoli/analysis/cohort/centralCarbonMetabolismScatter.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20180813.102120.814493__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 3

# panel C
python models/ecoli/analysis/single/massFractionSummary.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20180813.102120.814493__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 6 --generation 0 --seed 0

python models/ecoli/analysis/cohort/kinetics_flux_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20180813.102120.814493__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 6

python models/ecoli/analysis/cohort/centralCarbonMetabolismScatter.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20180813.102120.814493__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 6

# panels D, E and F
python models/ecoli/analysis/variant/metabolism_kinetic_objective_weight.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20180813.102120.814493__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/


## Supplemental figure 4
python models/ecoli/analysis/variant/meneSensitivity.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_F/20180814.113407.723416__SET_F_8_gens_8_seeds_9_menE_expression_values


## Supplemental figure 5
# made by JC from his experimental results
