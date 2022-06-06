# Script to run analysis to generate published figures in the operon paper
# Simulations in runscripts/paper/operon_paper/paper_runs.sh must be run to
# generate the data before performing analysis


## Figure 1
# Panel B
python models/ecoli/analysis/comparison/polycistronic_transcription.py \
out/20220602.133728__SET_E_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20220518.174543__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media


## Figure 2
# Panel B
python models/ecoli/analysis/comparison/doubling_time_histogram.py \
out/20220518.174606__SET_B_8_gens_128_seeds_operons_v1_with_glucose_minimal_media \
out/20220518.174543__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel C
python models/ecoli/analysis/comparison/mRNA_copy_numbers_growth_genes.py \
out/20220518.174606__SET_B_8_gens_128_seeds_operons_v1_with_glucose_minimal_media \
out/20220518.174543__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel E
# See Rend-seq repository

# Panel G
python models/ecoli/analysis/comparison/mRNA_copy_numbers_growth_genes.py \
out/20220526.173929__SET_C_8_gens_128_seeds_operons_v2_with_glucose_minimal_media \
out/20220518.174543__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel H
python models/ecoli/analysis/comparison/doubling_time_histogram.py \
out/20220526.173929__SET_C_8_gens_128_seeds_operons_v2_with_glucose_minimal_media \
out/20220518.174543__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media


## Figure 3
# Panel A (+Figure S2)
python models/ecoli/analysis/comparison/mRNA_copy_numbers.py \
out/20220526.173929__SET_C_8_gens_128_seeds_operons_v2_with_glucose_minimal_media \
out/20220518.174543__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel D (+Figure S3)
python models/ecoli/analysis/comparison/mRNA_copy_numbers.py \
out/20220518.174649__SET_D_8_gens_128_seeds_operons_v3_with_glucose_minimal_media \
out/20220518.174543__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel F, G
# See Rend-seq repository


## Figure 4
# Panel A
python models/ecoli/analysis/cohort/centralCarbonMetabolismScatter.py \
out/20220602.133728__SET_E_8_gens_128_seeds_operons_on_with_glucose_minimal_media

# Panel B
python models/ecoli/analysis/multigen/proteinCountsValidation.py \
out/20220602.133728__SET_E_8_gens_128_seeds_operons_on_with_glucose_minimal_media

# Panel C
python models/ecoli/analysis/comparison/mRNA_length_histogram.py \
out/20220602.133728__SET_E_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20220518.174543__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media


## Figure 5
# Panels A, B, C
python models/ecoli/analysis/comparison/coexpression_probabilities.py \
out/20220602.133728__SET_E_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20220518.174543__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panels E, F
python models/ecoli/analysis/comparison/excess_protein_monomers.py \
out/20220602.133728__SET_E_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20220518.174543__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

