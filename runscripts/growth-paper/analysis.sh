# Analysis of simulation results for each figure was performed with the commands below.
# Newly run simulations will have updated timestamps/paths.

# Figure 2a:
python runscripts/manual/analysisCohort.py out/20220308.182935__Long_lead_in_down_and_up_shifts_with_no_regulation/ -p growth_time_series
# Figure 2b:
python runscripts/manual/analysisCohort.py out/20220306.133934__Long_lead_in_down_and_up_shifts_with_regulation/ -p growth_time_series
# Figure 2c:
python runscripts/manual/analysisVariant.py out/20220121.200016__Amino_acid_combinations_in_media_without_regulation_or_charging/ -p doubling_time_histogram -o 6+gen-var0,3- --generation-path-range 6 25 --variant-path 0 3
# Figure 2d:
python runscripts/manual/analysisVariant.py out/20220116.130915__Amino_acid_combinations_in_media/ -p doubling_time_histogram -o 6+gen- --generation-path-range 6 25
# Figure 2e:
python runscripts/manual/analysisVariant.py out/20220116.130915__Amino_acid_combinations_in_media/ -p growth_trajectory  --generation-path-range 6 25
# Figure 2e:
python runscripts/manual/analysisVariant.py out/20220118.121155__Add_one_amino_acid_shift/ -p growth_trajectory  --generation-path-range 2 16
# Figure 2e:
python runscripts/manual/analysisVariant.py out/20220121.200028__Remove_one_amino_acid_shift/ -p growth_trajectory  --generation-path-range 4 8
# Figure 2e:
python runscripts/manual/analysisVariant.py out/20220123.144433__Conditions_without_regulation_or_charging/ -p growth_trajectory
# Figure 2e:
python runscripts/manual/analysisVariant.py out/20220123.144515__Conditions_with_regulation/ -p growth_trajectory
# Figure 2e:
runscripts/growth-paper/growth_rp_plot.py
# Figure 2f:
python runscripts/manual/analysisParca.py out/20220306.133934__Long_lead_in_down_and_up_shifts_with_regulation/ -p amino_acid_uptake_rates
# Figure 3a:
python runscripts/manual/analysisVariant.py out/20220117.063438__ppGpp_sensitivity -p growth_trajectory  --generation-path-range 2 8
# Figure 3a:
runscripts/growth-paper/growth_rp_plot.py
# Figure 3bcdefhi:
python runscripts/manual/analysisVariant.py out/20220117.063438__ppGpp_sensitivity -p ppgpp_conc -o 2+gen-minimal- --variant-path-range 0 10 --generation-path-range 2 8
# Figure 3gj:
python runscripts/manual/analysisVariant.py out/20220117.063438__ppGpp_sensitivity -p growth_trajectory -o 2+gen- --generation-path-range 2 8
# Figure 3g:
python runscripts/manual/analysisVariant.py out/20220117.215105__ppGpp_limitations_-_low_ppGpp/ -p growth_trajectory -o 2+gen- --generation-path-range 2 8
# Figure 3g:
runscripts/growth-paper/ppgpp_growth.py
# Figure 3j:
python runscripts/manual/analysisVariant.py out/20220120.060837__ppGpp_limitations_-_high_ppGpp/ -p growth_trajectory -o 2+gen- --generation-path-range 2 8
# Figure 3j:
python runscripts/manual/analysisVariant.py out/20220304.172940__ppGpp_limitations_with_ribosomes_at_high_ppGpp/ -p growth_trajectory -o 2+gen- --generation-path-range 2 8
# Figure 3j:
python runscripts/manual/analysisVariant.py out/20220304.172940__ppGpp_limitations_with_ribosomes_at_high_ppGpp,_no_ppGpp_translation_inhibition/ -p growth_trajectory -o 2+gen- --generation-path-range 2 8
# Figure 3j:
runscripts/growth-paper/ppgpp_growth.py
# Figure 4bef:
python runscripts/manual/analysisVariant.py out/20220119.081756__Remove_amino_acid_inhibition/ -p remove_aa_inhibition
# Figure 4b:
runscripts/growth-paper/remove-aa-inhib/conc_ki.py
# Figure 4c:
python runscripts/manual/analysisCohort.py out/20220315.074044__Remove_amino_acid_inhibition_large/ -p growth_time_series -o wt- -v0
# Figure 4c:
python runscripts/manual/analysisCohort.py out/20220315.074044__Remove_amino_acid_inhibition_large/ -p growth_time_series -o leuA- -v4
# Figure 4d:
python runscripts/manual/analysisVariant.py out/20220315.074044__Remove_amino_acid_inhibition_large/ -p aa_period
# Figure 4e:
runscripts/growth-paper/remove-aa-inhib/conc_ki.py
# Figure 4f:
runscripts/growth-paper/remove-aa-inhib/conc_ki.py
# Figure 5abc:
python runscripts/manual/analysisVariant.py out/20220116.130915__Amino_acid_combinations_in_media/ -p growth_trajectory  --generation-path-range 6 25
# Figure 5abc:
runscripts/growth-paper/growth_rp_plot.py
# Figure 5a:
python runscripts/manual/analysisVariant.py out/20220306.133934__Long_lead_in_down_and_up_shifts_with_regulation/ -p growth_trajectory
# Figure 5b:
python runscripts/manual/analysisVariant.py out/20220307.171936__Long_lead_in_down_and_up_shifts_without_mechanistic_translation_supply/ -p growth_trajectory
# Figure 5c:
python runscripts/manual/analysisVariant.py out/20220307.171901__Long_lead_in_down_and_up_shifts_without_ppGpp_regulation/ -p growth_trajectory
# Figure 5d:
python runscripts/manual/analysisCohort.py out/20220306.133934__Long_lead_in_down_and_up_shifts_with_regulation/ -p growth_time_series
# Figure 5e:
python runscripts/manual/analysisCohort.py out/20220307.171936__Long_lead_in_down_and_up_shifts_without_mechanistic_translation_supply/ -p growth_time_series
# Figure 5f:
python runscripts/manual/analysisCohort.py out/20220307.171901__Long_lead_in_down_and_up_shifts_without_ppGpp_regulation/ -p growth_time_series
# Figure 5ghi:
python runscripts/manual/analysisCohort.py out/20220307.171901__Long_lead_in_down_and_up_shifts_without_ppGpp_regulation/ -p growth_time_series -o no-ppgpp-
# Figure 5ghi:
python runscripts/manual/analysisCohort.py out/20220306.133934__Long_lead_in_down_and_up_shifts_with_regulation/ -p growth_time_series -o all-reg-
