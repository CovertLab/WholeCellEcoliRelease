# Scripts used to run analyses that generate the published figures in the operon
# paper
# You should run the simulations in runscripts/paper/operon_paper/paper_runs.sh
# prior to running these analyses. You will need to edit the timestamps in the
# directory names.


## Figure 1
# Panel B
python models/ecoli/analysis/comparison/polycistronic_transcription.py \
out/20220822.165041__SET_E_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20220822.164908__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media


## Figure 2
# Panel B
python models/ecoli/analysis/comparison/doubling_time_histogram.py \
out/20220822.164932__SET_B_8_gens_128_seeds_operons_v1_with_glucose_minimal_media \
out/20220822.164908__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel C
python models/ecoli/analysis/comparison/mRNA_copy_numbers_growth_genes.py \
out/20220822.164932__SET_B_8_gens_128_seeds_operons_v1_with_glucose_minimal_media \
out/20220822.164908__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel E
# See Rend-seq repository

# Panel G
python models/ecoli/analysis/comparison/mRNA_copy_numbers_growth_genes.py \
out/20220822.164954__SET_C_8_gens_128_seeds_operons_v2_with_glucose_minimal_media \
out/20220822.164908__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel H
python models/ecoli/analysis/comparison/doubling_time_histogram.py \
out/20220822.164954__SET_C_8_gens_128_seeds_operons_v2_with_glucose_minimal_media \
out/20220822.164908__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media


## Figure 3
# Panel A
python models/ecoli/analysis/comparison/mRNA_copy_numbers.py \
out/20220822.164954__SET_C_8_gens_128_seeds_operons_v2_with_glucose_minimal_media \
out/20220822.164908__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel E
python models/ecoli/analysis/parca/corrected_rnaseq_read_counts.py \
out/20220822.165018__SET_D_8_gens_128_seeds_operons_v3_with_glucose_minimal_media

# Panel F (+ Figure S2A)
python models/ecoli/analysis/comparison/mRNA_copy_numbers.py \
out/20220822.165018__SET_D_8_gens_128_seeds_operons_v3_with_glucose_minimal_media \
out/20220822.164908__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panels H, I
# See Rend-seq repository


## Figure 4
# Panel A
python models/ecoli/analysis/comparison/mRNA_length_histogram.py \
out/20220822.165041__SET_E_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20220822.164908__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel B
python models/ecoli/analysis/comparison/mRNA_counts_histogram.py \
out/20220822.165041__SET_E_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20220822.164908__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel C
python models/ecoli/analysis/comparison/mRNA_mass_histogram.py \
out/20220822.165041__SET_E_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20220822.164908__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media


## Figure 5
# Panel A (+ Table S3)
python models/ecoli/analysis/comparison/coexpression_probabilities.py \
out/20220822.165041__SET_E_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20220822.164908__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panels B, C, D (+ Table S4)
python models/ecoli/analysis/comparison/protein_stoichiometry.py \
out/20220822.165041__SET_E_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20220822.164908__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panels E, F, G (+ Table S5)
python models/ecoli/analysis/comparison/excess_protein_monomers.py \
out/20220822.165041__SET_E_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20220822.164908__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media


## Figure S2
# Panel A
python models/ecoli/analysis/comparison/mRNA_copy_numbers.py \
out/20220822.165018__SET_D_8_gens_128_seeds_operons_v3_with_glucose_minimal_media \
out/20220822.164908__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel B
# See Rend-seq repository


## Figure S3
# Panel A
python models/ecoli/analysis/comparison/doubling_time_histogram.py \
out/20220822.165041__SET_E_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20220822.164908__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panels B, C
python models/ecoli/analysis/comparison/proteomics_fluxomics_comparison.py \
out/20220822.165041__SET_E_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20220822.164908__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media


## Table S1
python models/ecoli/analysis/parca/mRNA_transcript_table.py \
out/20220822.165041__SET_E_8_gens_128_seeds_operons_on_with_glucose_minimal_media
