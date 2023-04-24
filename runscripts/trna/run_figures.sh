# Scripts used to run analyses that generate the figures in paper:
# "Whole-cell modeling of E. coli confirms that in vitro tRNA
# aminoacylation measurements are insufficient to support cell growth
# and predicts a positive feedback mechanism regulating arginine
# biosynthesis".

# The simulations described in runscripts/paper/trna/run_simulations.sh
# must be run prior to running the analyses in this file.


## Figure 2
# Panel A
python models/ecoli/analysis/variant/trna_synthetase_kinetics_variant_doubling.py sims

# Panel B
python models/ecoli/analysis/cohort/mass_cell_cycle.py -v 2 sims

# Panel C
python models/ecoli/analysis/cohort/mass_cell_cycle.py -v 2 sims

# Panel E
python models/ecoli/analysis/cohort/trna_synthetase_average.py -v 1 sims

# Panel F
python models/ecoli/analysis/cohort/trna_synthetase_distribution.py -v 1 sims


## Figure 3
# Panel B
python models/ecoli/analysis/variant/trna_synthetase_kinetics_variant_doubling.py sims

# Panel C
python models/ecoli/analysis/variant/trna_charging_views.py sims

# Panel D
python models/ecoli/analysis/variant/trna_charging_validation.py sims

# Panel E
python models/ecoli/analysis/variant/trna_charging_validation.py sims

# Panel F
python models/ecoli/analysis/variant/trna_charging_validation.py sims

# Panel G
python models/ecoli/analysis/variant/trna_charging_validation.py sims

# Panel H
python models/ecoli/analysis/parca/trna_synthetase_kcats.py sims
    

## Figure 4
# Panel A
python models/ecoli/analysis/variant/doubling_time_beeswarm.py sims

# Panel B
python models/ecoli/analysis/variant/HisRS_kcat_impact.py sims

# Panel C
python models/ecoli/analysis/variant/HisRS_kcat_impact.py sims

# Panel D
python models/ecoli/analysis/variant/HisRS_kcat_impact.py sims

# Panel E
python models/ecoli/analysis/variant/HisRS_kcat_impact.py sims

# Panel F
python models/ecoli/analysis/variant/synthetase_ribosome_relation.py sims

# Panel H
python models/ecoli/analysis/variant/synthetase_ribosome_relation_sweep.py sims


# Figure 5
# Panel A
python models/ecoli/analysis/variant/doubling_time_bars.py sims

# Panel B
python models/ecoli/analysis/variant/mass_single_lineage.py sims

# Panel C
python models/ecoli/analysis/variant/mass_cell_cycle_fold_change.py sims

# Panel D
python models/ecoli/analysis/cohort/trna_charged_fractions.py -v 6 sims

# Panel E
python models/ecoli/analysis/variant/argRS_and_substrates.py sims

# Panel F
python models/ecoli/analysis/variant/arginine_biosynthesis.py sims

# Panel G
python models/ecoli/analysis/variant/argA_counts.py sims

# Panel H
python models/ecoli/analysis/variant/argA_ribosome_profile.py sims

# Panel I
python models/ecoli/analysis/variant/argA_ribosome_terminations.py sims


## Figure 6
# Panel A
python models/ecoli/analysis/variant/argRS_kcat_impact.py sims


## Figure S2
python models/ecoli/analysis/cohort/trna_synthetase_distribution.py -v 1 sims


## Figure S3
python models/ecoli/analysis/variant/trna_synthetase_dynamic_range_sweep.py out


## Figure S4
# Panel A
python models/ecoli/analysis/cohort/trna_forms.py -v 0 sims

# Panel B
python models/ecoli/analysis/cohort/trna_charged_fractions_validation.py -v 0 sims


## Figure S5
# Panel A
python models/ecoli/analysis/cohort/premature_ribosome_termination.py -v 6 sims

# Panel B
python models/ecoli/analysis/cohort/premature_ribosome_termination.py -v 6 sims

# Panel C
python models/ecoli/analysis/variant/rare_arginine_codons_expression.py sims

# Panel D
python models/ecoli/analysis/variant/rare_arginine_codons_expression.py sims

# Panel E
python models/ecoli/analysis/variant/argA_codon_experiment.py sims
