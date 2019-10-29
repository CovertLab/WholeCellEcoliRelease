# Script to run analysis to generate published figures
# When reproducing, /scratch/PI/mcovert/wc_ecoli/paper will be the out directory
# Simulations in runscripts/paper/paper_runs.sh must be run to generate the data before performing analysis


## Figure 1
# Bottom Simulations Left to Right, Top to bottom:
python models/ecoli/analysis/multigen/growth_dynamics_grid.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_B/20190908.151927__SET_B_9_gens_8_seeds_shift_to_plus_AA_without_growth_noise_with_D_period/ \
--seed 0


## Figure 2
# Panel A
python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_L/20190910.155143__SET_L_4_gens_256_seeds_3_conditions_unfit_ribosome_and_rna_poly_expression/

# Panel B
python models/ecoli/analysis/variant/param_sensitivity.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_K/SET_K_sensitivity/

# Panel C
python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_D/20190908.203751__SET_D2_4_gens_256_seeds,_unfit_ribosome_expression/

python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_D/20190909.171540__SET_D3_4_gens_256_seeds,_unfit_rna_poly_expression/

python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_D/20190909.224043__SET_D4_4_gens_256_seeds,_unfit_ribosome_and_rna_poly_expression,_adjusted_RNase_expression/

# Panel D
python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20190906.151026__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/

# Panel E
# Simulation values were obtained from pausing
# reconstruction/ecoli/fit_sim_data_1.py after line 381, reading abundances of
# RNA Polymerase and Ribosome subunits from the bulkContainer object, and
# calculating the number of full complexes that can be formed from their stoichiometries.

# Panel F
python paper/figure_s2_investigation/mrna_expression/mrna_expression.py

# Panel G
python models/ecoli/analysis/cohort/rnap_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.090047__SET_J03_alternate_RNA_seq

python models/ecoli/analysis/cohort/ribosome_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.090047__SET_J03_alternate_RNA_seq

# Panel H
python models/ecoli/analysis/variant/growth_condition_comparison_validation_parallel.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20190906.151026__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/

## Figure 3
# panel A
python models/ecoli/analysis/cohort/limited_metabolites.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_A/20190910.135308__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period

python models/ecoli/analysis/cohort/limited_metabolites.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_G/20190919.184128__SET_G_32_gens_8_seeds_basal_without_parameter_adjustments

# panels B and G
python models/ecoli/analysis/variant/glc_yield.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_I/20190910.155355__SET_I_kinetic_constraint_factorial_design/

# panel C
python models/ecoli/analysis/cohort/centralCarbonMetabolismScatter.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_M/20190919.104605__SET_M1_old_kinetic_constraints

# panel D
python models/ecoli/analysis/variant/flux_sensitivity.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_H/20190906.115628__SET_H_1_gen_flux_sensitivity/

# panel E
python models/ecoli/analysis/variant/kinetic_objective_interactions.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_I/20190910.155355__SET_I_kinetic_constraint_factorial_design/

# panel F
python models/ecoli/analysis/variant/kinetic_objective_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_I/20190910.155355__SET_I_kinetic_constraint_factorial_design/

# panel H
python models/ecoli/analysis/cohort/centralCarbonMetabolismScatter.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_M/20190919.162726__SET_M2_new_kinetic_constraints

# panel I
python models/ecoli/analysis/single/massFractionSummary.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_I/20190910.155355__SET_I_kinetic_constraint_factorial_design/ \
--variant 47 --seed 0 --gen 0

# panel J
python models/ecoli/analysis/cohort/kinetics_flux_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_M/20190919.162726__SET_M2_new_kinetic_constraints


## Figure 4
# panel A
python models/ecoli/analysis/multigen/proteinCountsValidation.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_A/20190910.135308__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/ \
--seed 0

# panel B
python models/ecoli/analysis/cohort/expression_dynamics.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_A/20190910.135308__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/

# panel C and D
python models/ecoli/analysis/multigen/subgenerational_transcription.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_A/20190910.135308__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/ \
--seed 0

# panel E
python models/ecoli/analysis/multigen/functionalUnits.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_A/20190910.135308__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/ \
--seed 0

# panel F
python models/ecoli/analysis/multigen/pabx_limitations.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_A/20190910.135308__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/ \
--seed 5


## Figure 5
# panel A
python paper/steadyStateAnalysis/main.py

# panel B - experimental results

# panel C - experimental results


## Supplemental figure 1
# panel A
python models/ecoli/analysis/multigen/growth_dynamics_panel.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_B/20190908.151927__SET_B_9_gens_8_seeds_shift_to_plus_AA_without_growth_noise_with_D_period/ \
--seed 0

# panels B and C
python models/ecoli/analysis/multigen/environmental_shift_fluxes.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_B/20190908.151927__SET_B_9_gens_8_seeds_shift_to_plus_AA_without_growth_noise_with_D_period/ \
--seed 0


## Supplemental figure 2
# Panel A
# Input Comparison
python paper/figure_s2_investigation/mrna_expression/mrna_expression.py
python paper/figure_s2_investigation/translation_efficiency/translation_efficiency.py
python paper/figure_s2_investigation/mrna_half_lives/mrna_half_lives.py

# RNA Polymerase (ordered from column 0 to 11)
python models/ecoli/analysis/cohort/rnap_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_L/20190910.155143__SET_L_4_gens_256_seeds_3_conditions_unfit_ribosome_and_rna_poly_expression

python models/ecoli/analysis/cohort/rnap_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190923.193545__SET_J01_alternate_RNA_mass_per_cell

python models/ecoli/analysis/cohort/rnap_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190923.223401__SET_J02_alternate_mRNA_mass_per_cell

python models/ecoli/analysis/cohort/rnap_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.090047__SET_J03_alternate_RNA_seq

python models/ecoli/analysis/cohort/rnap_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.111334__SET_J04_alternate_translation_efficiency

python models/ecoli/analysis/cohort/rnap_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.134144__SET_J05_alternate_R-protein_half_lives

python models/ecoli/analysis/cohort/rnap_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.180126__SET_J06_alternate_protein_mass_per_cell

python models/ecoli/analysis/cohort/rnap_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.200659__SET_J07_variable_transcription_elongation_rate

python models/ecoli/analysis/cohort/rnap_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.220818__SET_J08_max_RNA_polymerase_activity

python models/ecoli/analysis/cohort/rnap_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.234820__SET_J09_mRNA_half_life_fitting

python models/ecoli/analysis/cohort/rnap_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190925.014132__SET_J10_variable_translation_elongation_rate

python models/ecoli/analysis/cohort/rnap_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190925.070638__SET_J11_alternate_ribosome_activity

# Ribosome (ordered from column 0 to 11)
python models/ecoli/analysis/cohort/ribosome_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_L/20190910.155143__SET_L_4_gens_256_seeds_3_conditions_unfit_ribosome_and_rna_poly_expression

python models/ecoli/analysis/cohort/ribosome_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190923.193545__SET_J01_alternate_RNA_mass_per_cell

python models/ecoli/analysis/cohort/ribosome_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190923.223401__SET_J02_alternate_mRNA_mass_per_cell

python models/ecoli/analysis/cohort/ribosome_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.090047__SET_J03_alternate_RNA_seq

python models/ecoli/analysis/cohort/ribosome_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.111334__SET_J04_alternate_translation_efficiency

python models/ecoli/analysis/cohort/ribosome_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.134144__SET_J05_alternate_R-protein_half_lives

python models/ecoli/analysis/cohort/ribosome_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.180126__SET_J06_alternate_protein_mass_per_cell

python models/ecoli/analysis/cohort/ribosome_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.200659__SET_J07_variable_transcription_elongation_rate

python models/ecoli/analysis/cohort/ribosome_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.220818__SET_J08_max_RNA_polymerase_activity

python models/ecoli/analysis/cohort/ribosome_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.234820__SET_J09_mRNA_half_life_fitting

python models/ecoli/analysis/cohort/ribosome_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190925.014132__SET_J10_variable_translation_elongation_rate

python models/ecoli/analysis/cohort/ribosome_validation.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190925.070638__SET_J11_alternate_ribosome_activity

# Doubling Time (ordered from column 0 to 11)
python models/ecoli/analysis/cohort/doubling_times_violin.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_L/20190910.155143__SET_L_4_gens_256_seeds_3_conditions_unfit_ribosome_and_rna_poly_expression/

python models/ecoli/analysis/cohort/doubling_times_violin.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190923.193545__SET_J01_alternate_RNA_mass_per_cell

python models/ecoli/analysis/cohort/doubling_times_violin.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190923.223401__SET_J02_alternate_mRNA_mass_per_cell

python models/ecoli/analysis/cohort/doubling_times_violin.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.090047__SET_J03_alternate_RNA_seq

python models/ecoli/analysis/cohort/doubling_times_violin.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.111334__SET_J04_alternate_translation_efficiency

python models/ecoli/analysis/cohort/doubling_times_violin.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.134144__SET_J05_alternate_R-protein_half_lives

python models/ecoli/analysis/cohort/doubling_times_violin.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.180126__SET_J06_alternate_protein_mass_per_cell

python models/ecoli/analysis/cohort/doubling_times_violin.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.200659__SET_J07_variable_transcription_elongation_rate

python models/ecoli/analysis/cohort/doubling_times_violin.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.220818__SET_J08_max_RNA_polymerase_activity

python models/ecoli/analysis/cohort/doubling_times_violin.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190924.234820__SET_J09_mRNA_half_life_fitting

python models/ecoli/analysis/cohort/doubling_times_violin.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190925.014132__SET_J10_variable_translation_elongation_rate

python models/ecoli/analysis/cohort/doubling_times_violin.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_J/20190925.070638__SET_J11_alternate_ribosome_activity

# Panel C
python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20190906.151026__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/ \
--variant_index 1

python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20190906.151026__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/ \
--variant_index 2

# panel D
python models/ecoli/analysis/variant/rna_protein_ratio_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20190906.151026__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/

python models/ecoli/analysis/variant/rna_protein_ratio_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_L/20190910.155143__SET_L_4_gens_256_seeds_3_conditions_unfit_ribosome_and_rna_poly_expression/

# panel E
python models/ecoli/analysis/variant/adder_sizer.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20190906.151026__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/

# panel F
python models/ecoli/analysis/variant/adder_sizer_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20190906.151026__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/


## Supplemental figure 3
# panel A
python models/ecoli/analysis/single/massFractionSummary.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20190917.224926__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 0 --generation 0 --seed 0

# panel B
python models/ecoli/analysis/single/massFractionSummary.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20190917.224926__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 6 --generation 0 --seed 0

# panel C
python models/ecoli/analysis/cohort/kinetics_flux_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20190917.224926__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 0

# panel D
python models/ecoli/analysis/cohort/kinetics_flux_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20190917.224926__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 6

# panels E and F
python models/ecoli/analysis/variant/metabolism_kinetic_objective_weight.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20190917.224926__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/

# panel G - experimental results

# panel H - same as figure 3D
python models/ecoli/analysis/variant/flux_sensitivity.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_H/20190906.115628__SET_H_1_gen_flux_sensitivity/



## Supplemental figure 4
python models/ecoli/analysis/variant/subgen_expression.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_F/20190918.093118__SET_F_16_gens_8_seeds_9_pabB_expression_values/


## Supplemental figure 5 - experimental results
