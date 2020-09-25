'''Generate development analysis plots for colony simulations

For usage information, execute:

	python analysis.py -h
'''


from vivarium_cell.analysis.analyze import Analyzer


TIMESERIES_CONFIG = {
	'skip_paths': [
		('boundary', 'wcecoli_fields_null'),
	],
}
SNAPSHOTS_CONFIG = {
	'include_fields': ['nitrocefin', 'GLC'],
	'field_label_size': 54,
	'default_font_size': 54,
}
TAGS_CONFIG = {
	'tag_label_size': 54,
	'default_font_size': 54,
}


def main():
	analyzer = Analyzer(
		timeseries_config=TIMESERIES_CONFIG,
		snapshots_config=SNAPSHOTS_CONFIG,
	)
	analyzer.run()


if __name__ == '__main__':
	main()
