'''Generate development analysis plots for colony simulations

For usage information, execute:

    python analysis.py -h
'''

from vivarium.analysis.analyze import Analyzer


TIMESERIES_CONFIG = {
	'remove_zeros': True,
	'skip_paths': [
		('boundary', 'wcecoli_fields_null'),
	],
}
SNAPSHOTS_CONFIG = {
	'include_fields': ['nitrocefin', 'GLC'],
}


def main():
	analyzer = Analyzer(
		timeseries_config=TIMESERIES_CONFIG,
		snapshots_config=SNAPSHOTS_CONFIG,
	)
	analyzer.run()


if __name__ == '__main__':
    main()
