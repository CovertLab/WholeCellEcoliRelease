from __future__ import absolute_import, division, print_function

# Enable segmentation and other fault handling for tracebacks
# noinspection PyCompatibility
import faulthandler; faulthandler.enable()

# Set OPENBLAS_NUM_THREADS to 1 to avoid threading bugs in OpenBLAS
# that can lead to irreproducible results across different systems. It may have
# performance implications and could be removed if a library test is created to
# make sure the installed version of OpenBLAS is not dependent on the number of
# threads allowed. Must be set before numpy (and maybe scipy) is loaded.
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

from .initRawData import InitRawDataTask
from .fitSimData import FitSimDataTask
from .variantSimData import VariantSimDataTask
from .simulation import SimulationTask
from .initRawValidationData import InitRawValidationDataTask
from .initValidationData import InitValidationDataTask
from .simulationDaughter import SimulationDaughterTask
from .analysisSingle import AnalysisSingleTask
from .analysisMultiGen import AnalysisMultiGenTask
from .analysisCohort import AnalysisCohortTask
from .analysisVariant import AnalysisVariantTask
from .analysisParca import AnalysisParcaTask
from .buildCausalityNetwork import BuildCausalityNetworkTask
from .parca import ParcaTask
from .writeJson import WriteJsonTask
