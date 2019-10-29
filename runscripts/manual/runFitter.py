"""
Run the Fitter. The output goes into the named subdirectory of wcEcoli/out/,
defaulting to "manual".

TODO: Share lots of code with fw_queue.py and AnalysisPaths.py.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import
from __future__ import division

import os

from wholecell.fireworks.firetasks import FitSimDataTask
from wholecell.fireworks.firetasks import InitRawDataTask
from wholecell.fireworks.firetasks import InitRawValidationDataTask
from wholecell.fireworks.firetasks import InitValidationDataTask
from wholecell.fireworks.firetasks import SymlinkTask
from wholecell.utils import constants, scriptBase
from wholecell.utils import filepath as fp


class RunFitter(scriptBase.ScriptBase):
	"""Runs the Fitter, aka simulation parameter calculator."""

	def define_parameters(self, parser):
		super(RunFitter, self).define_parameters(parser)

		# NOTE: Don't name this arg sim_dir since that makes parse_args() look
		# for an existing sim_dir directory.
		parser.add_argument('sim_outdir', nargs='?', default='manual',
			help='The simulation "out/" subdirectory to write to.'
				 ' Default = "manual".')
		parser.add_argument('--timestamp', action='store_true',
			help='Timestamp the given `sim_outdir`, transforming e.g.'
				 ' "Fast run" to "20190514.135600__Fast_run".')
		parser.add_argument('-c', '--cpus', type=int, default=1,
			help='The number of CPU processes to use. Default = 1.')
		parser.add_argument('--cached', action='store_true',
			help='Copy a cached "' + constants.SERIALIZED_FIT1_FILENAME
				 + '" file instead of generating it.')
		parser.add_argument('-d', '--debug', action='store_true',
			help="Enable Fitter debugging: Fit only one arbitrarily-chosen"
				 " transcription factor for a faster debug cycle. Don't use it"
				 " for an actual simulation.")
		parser.add_argument(
			'--disable-ribosome-fitting',
			action='store_true',
			help= "If set, ribosome expression will not be fit to protein synthesis demands."
			)
		parser.add_argument(
			'--disable-rnapoly-fitting',
			action='store_true',
			help= "If set, RNA polymerase expression will not be fit to protein synthesis demands."
			)
		parser.add_argument(
			'--variable-elongation-transcription',
			action='store_true',
			help= "If set, one elongation rate will be used for transcription."
			)
		parser.add_argument(
			'--variable-elongation-translation',
			action='store_true',
			help= "If set, one elongation rate will be used for translation"
			)
		parser.add_argument(
			'--no-expression-adjustment',
			action='store_false',
			help= "If set, some RNA and protein expression parameters will not be adjusted."
			)
		parser.add_argument(
			'--adjust-rnase-expression',
			action='store_true',
			help= "If set, adjusts the expression of all RNase mRNA lower."
			)
		parser.add_argument(
			'--alternate-mass-fraction-mrna',
			action='store_true',
			help="If set, allocates smaller mass fraction for mRNA."
		)
		parser.add_argument(
			'--alternate-mass-fraction-protein',
			action='store_true',
			help="If set, allocates larger mass fraction for protein."
			)
		parser.add_argument(
			'--alternate-mass-fraction-rna',
			action='store_true',
			help="If set, allocates smaller mass fraction for RNA."
		)
		parser.add_argument(
			'--alternate-r-protein-degradation',
			action='store_true',
			help= "If set, r-protein degradation will set to fast value."
			)
		parser.add_argument(
			'--alternate-ribosome-activity',
			action='store_true',
			help="Alternate ribosome active fraction: 85 percent active. Default = 80 percent active."
		)
		parser.add_argument(
			'--alternate-rna-half-life',
			action='store_true',
			help="If set, alternate RNA half life input will be used."
		)
		parser.add_argument(
			'--alternate-rna-seq',
			action='store_true',
			help="If set, alternate RNA-seq (Covert 2004) input will be used."
			)
		parser.add_argument(
			'--alternate-translation-efficiency',
			action='store_true',
			help="Alternate translation efficiency described by Mohammad et"
				 "al. 2019. Default = described by Li et al. 2014."
		)
		parser.add_argument(
			'--disable-measured-protein-deg',
			action='store_true',
			help="If set, does not use any measured protein degradation rates"
			     " and defaults to the N-end rule."
		)
		parser.add_argument(
			'--disable-ribosome-activity-fix',
			action='store_true',
			help="If set, disables ribosome activity fix."
		)
		parser.add_argument(
			'--disable-rnap-fraction-increase',
			action='store_true',
			help="If set, disables doubling-time-dependent RNAP fraction increase."
		)
		parser.add_argument(
			'--max-rnap-activity',
			action='store_true',
			help="If set, RNA polymerase activity will be set to 100 percent."
		)
		parser.add_argument(
			'--mrna-half-life-fitting',
			action='store_true',
			help="If set, mRNA half lives will be fit to transcription demands."
		)
		parser.add_argument(
			'--rnapoly-activity-fitting',
			action='store_true',
			help="If set, RNA polymerase activity will be fit to transcription demands."
		)
		parser.add_argument(
			'--save-cell-specs',
			action='store_true',
			help="If set, saves cell specs."
		)
		parser.add_argument(
			'--write-translation-efficiencies',
			action='store_true',
			help="If set, writes out translation efficiencies."
		)

	def parse_args(self):
		args = super(RunFitter, self).parse_args()

		if args.timestamp:
			args.sim_outdir = fp.timestamp() + '__' + args.sim_outdir.replace(
				' ', '_')

		args.sim_path = fp.makedirs(fp.ROOT_PATH, "out", args.sim_outdir)
		return args

	def run(self, args):
		kb_directory = fp.makedirs(args.sim_path, "kb")
		raw_data_file = os.path.join(kb_directory, constants.SERIALIZED_RAW_DATA)
		sim_data_file = os.path.join(kb_directory, constants.SERIALIZED_FIT1_FILENAME)
		cell_specs_file = os.path.join(kb_directory, constants.SERIALIZED_CELL_SPECS)
		cached_sim_data_file = os.path.join(
			fp.ROOT_PATH, 'cached', constants.SERIALIZED_FIT1_FILENAME)
		most_fit_filename = os.path.join(
			kb_directory, constants.SERIALIZED_SIM_DATA_MOST_FIT_FILENAME)
		raw_validation_data_file = os.path.join(
			kb_directory, constants.SERIALIZED_RAW_VALIDATION_DATA)
		validation_data_file = os.path.join(
			kb_directory, constants.SERIALIZED_VALIDATION_DATA)

		if args.debug or args.cached:
			print "{}{}Fitter".format(
				'DEBUG ' if args.debug else '',
				'CACHED ' if args.cached else '',
				)

		tasks = [
			InitRawDataTask(
				output=raw_data_file,
				),

			FitSimDataTask(
				fit_level=1,
				input_data=raw_data_file,
				output_data=sim_data_file,
				cached=args.cached,  # bool
				cached_data=cached_sim_data_file,  # cached file to copy
				cpus=args.cpus,
				debug=args.debug,
				disable_ribosome_capacity_fitting=args.disable_ribosome_fitting,
				disable_rnapoly_capacity_fitting=args.disable_rnapoly_fitting,
				variable_elongation_transcription=args.variable_elongation_transcription,
				variable_elongation_translation=args.variable_elongation_translation,
				rnapoly_activity_fitting=args.rnapoly_activity_fitting,
				mrna_half_life_fitting=args.mrna_half_life_fitting,
				max_rnap_activity=args.max_rnap_activity,
				adjust_rna_and_protein_parameters=args.no_expression_adjustment,
				adjust_rnase_expression=args.adjust_rnase_expression,
				disable_measured_protein_deg=args.disable_measured_protein_deg,
				alternate_mass_fraction_protein=args.alternate_mass_fraction_protein,
				alternate_mass_fraction_rna=args.alternate_mass_fraction_rna,
				alternate_mass_fraction_mrna=args.alternate_mass_fraction_mrna,
				alternate_r_protein_degradation=args.alternate_r_protein_degradation,
				alternate_rna_seq=args.alternate_rna_seq,
				alternate_rna_half_life=args.alternate_rna_half_life,
				alternate_translation_efficiency=args.alternate_translation_efficiency,
				alternate_ribosome_activity=args.alternate_ribosome_activity,
				disable_rnap_fraction_increase=args.disable_rnap_fraction_increase,
				disable_ribosome_activity_fix=args.disable_ribosome_activity_fix,
				save_cell_specs=args.save_cell_specs,
				cell_specs_file=cell_specs_file,
				write_translation_efficiencies=args.write_translation_efficiencies
				),

			SymlinkTask(
				to=constants.SERIALIZED_FIT1_FILENAME,
				link=most_fit_filename,
				overwrite_if_exists=True,
				),

			InitRawValidationDataTask(
				output=raw_validation_data_file,
				),

			InitValidationDataTask(
				validation_data_input=raw_validation_data_file,
				knowledge_base_raw=raw_data_file,
				output_data=validation_data_file,
				),
			]
		for task in tasks:
			task.run_task({})

		print '\n\t'.join(['Wrote', raw_data_file, sim_data_file,
			most_fit_filename, raw_validation_data_file, validation_data_file])


if __name__ == '__main__':
	script = RunFitter()
	script.cli()
