# Running the model

See the top level [README](../README.md) for general instructions and [docs/README](README.md) for details on setup and running
and tradeoffs between different ways to run the model.


## Fireworks

See [Setting up to run FireWorks](wholecell/fireworks/README.md) for setup to run a FireWorks workflow of cell simulations and analysis plots.


## Using the Manual Runscripts

[TODO] Rewrite this. Meanwhile, see the top level [README](../README.md).

You can accomplish largely the same thing that running via Fireworks does by running:

1. `runParca.py` -- creates the `kb/*` sim data files for the control case
2. `runSim.py` -- creates the modified modified sim data files for requested variants and runs the sims
3. `analysisCohort.py`, `analysisMultigen.py`, `analysisSingle.py`, `analysisVariant.py` -- runs the plots

Except:

* It runs directly without FireWorks or MongoDB.
* You can run it in a debugger.
* You can run just the parts you want and rerun them as needed _if_ you attend to the overall input and output files. The manual scripts don't automate dependency management. It's on you to rerun code if things change, runSim before analysis, or delete runSim output before running it again. (That last part should be improved! Also note that some analysis scripts get confused if the sim runs are more varied than expected. See [https://github.com/CovertLab/wcEcoli/issues/199](#199).)
* The scripts have command line help and are more rigorous at parsing args than fw_queue is at reading environment variables. Argparse provides helpful features like the ability to use any unambiguous option name prefix so you can use `--var` in place of `--variant_index`.
* The scripts do smart defaulting. runParca defaults to creating the sim dir `out/manual/` rather than a timestamped directory name. ßThe others default to finding the most likely sim dir, i.e. the latest timestamped subdirectory of `out/` or else the alphabetically first subdirectory of `out/`. You can pass in a directory name that's absolute, or relative to `wcEcoli/`, or relative to `wcEcoli/out/`.
* Analysis scripts can accept a parameter like `--variant_index 1` to look for a subdirectory like `condition_000001/`.
* You can run a single analysis script via

       python models/ecoli/analysis/cohort/growthDynamics.py

   or run one or more of them like this:

       runscripts/manual/analysisCohort.py growthDynamics transcriptFrequency
* The manual scripts don't compress any output files.
