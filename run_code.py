# import os
# import os.path
# import cPickle
# import wholecell.kb.KnowledgeBase

# # Create output directory
# outDir = "data/fixtures"
# if not os.path.exists(outDir):
#     os.makedirs(outDir)

# # Construct KB fixture
# kb = wholecell.kb.KnowledgeBase.KnowledgeBase(dataFileDir = "data/parsed", seqFileName = "data/raw/sequence.txt")
# cPickle.dump(kb, open(os.path.join(outDir, "KnowledgeBase.cPickle"), "w"), protocol = cPickle.HIGHEST_PROTOCOL)

import cPickle
import os
import os.path
import wholecell.kb.KnowledgeBase
import wholecell.sim.Simulation
import wholecell.sim.logger.Shell
import wholecell.sim.logger.DiskPyTables

# Load KB fixture
kbDir = "data/fixtures"
kb = cPickle.load(open(os.path.join(kbDir, "KnowledgeBase.cPickle"), "r"))

# Initialize simulation
sim = wholecell.sim.Simulation.Simulation(kb)
sim.setOptions({"seed": 10, "lengthSec": 100})

# Fit simulation
import wholecell.util.Fitter as Fitter
Fitter.Fitter.FitSimulation(sim, kb)
sim.calcInitialConditions()

# Initialize loggers
shellLogger = wholecell.sim.logger.Shell.Shell()

sim.run([shellLogger])
