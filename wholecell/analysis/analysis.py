
# TODO: port analysis scripts to this class
# TODO: simulation class method that calls analysis
# TODO: defaults for figure sizes, formatting, ...?

class AnalysisException(Exception):
	pass

# TODO: determine how plotOutDir will interact with cohort jobs

class Analysis(object):
	_interactiveOnly = False

	def __init__(self, simOutDir, plotOutDir = None, interactiveMode = False):
		if interactiveMode is False and self._interactiveOnly:
			raise AnalysisException(
				"{} must be ran in interactive mode.".format(self.__name__)
				)

		if plotOutDir is None and interactiveMode is False:
			raise AnalysisException("Analysis scripts not ran in interactive" +
				" mode must assign an output directory.")
		
		self._simOutDir = simOutDir
		self._plotOutDir = plotOutDir
		self._interactiveMode = interactiveMode

		self._filenamesUsed = []


	def run(self):
		raise NotImplementedError()


	def _loadKB(self, kbVersion = None):
		# load most fit KB by default
		raise NotImplementedError()


	def _loadDescription(self):
		raise NotImplementedError()


	def _loadGitDiff(self):
		raise NotImplementedError()


	def _loadGitHash(self):
		raise NotImplementedError()


	def _saveFigure(self, figure, subtitle = None):
		# append to self._filenamesUsed and raise if a filename is used twice
		raise NotImplementedError()


	@classmethod
	def interactiveOnly(cls):
		return cls._interactiveOnly


	def figureFilenames(self):
		return self._filenamesUsed


class AnalysisSingle(Analysis):
	def _outputDirectory(self):
		raise NotImplementedError()


class AnalysisCohort(Analysis):
	def _outputDirectories(self):
		raise NotImplementedError()
