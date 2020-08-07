# Running the model

See the top level [README](../README.md) for general instructions and [README](README.md) for details on setup and running.


## FireWorks

See [Setting up to run FireWorks](wholecell/fireworks/README.md) for setup to run a FireWorks workflow of cell simulations and analysis plots.


## Using the Manual Runscripts

See the top level [README](../README.md) for more details.

You can accomplish largely the same thing as a FireWorks run via:

1. `runFitter.py` -- creates the `kb/*` sim data files for the control case
2. `runSim.py` -- creates the modified modified sim data files for requested variants and runs the sims
3. `analysisCohort.py`, `analysisMultigen.py`, `analysisSingle.py`, `analysisVariant.py` -- runs the plots

In comparison:

* It runs directly without FireWorks or MongoDB.
* You can run it in a debugger.
* You can run just the parts you want and rerun them as needed _if_ you attend to the overall input and output files. The manual scripts don't automate dependency management. It's on you to rerun code if things change, runSim before analysis, or delete runSim output before running it again. (That last part should be improved! Also note that some analysis scripts get confused if the sim runs are more varied than expected. See [https://github.com/CovertLab/wcEcoli/issues/199](#199).)
* The scripts have command line help and are more rigorous at parsing args than fw_queue is at reading environment variables. Argparse provides helpful features like the ability to use any unambiguous option name prefix so you can use `--var` in place of `--variant_index`.
* The scripts do smart defaulting. runFitter defaults to creating the sim dir `out/manual/` rather than a timestamped directory name. (Would anyone prefer a timestamped directory?) The others default to finding the most likely sim dir, i.e. the latest timestamped subdirectory of `out/` or else the alphabetically first subdirectory of `out/`. You can pass in a directory name that's absolute, or relative to `wcEcoli/`, or relative to `wcEcoli/out/`.
* runSim supports variants but not multiple initial simulations, multiple generations, or multiple daughters per generation. These features are not hard to add if needed.
* The only parallelism supported is runFitter's `-cpus` arg and the `WC_ANALYZE_FAST` environment variable to the analysis scripts (which ought to be replaced with a command line arg).
* Analysis scripts can accept a parameter like `--variant_index 1` to look for a subdirectory like `condition_000001/`. (Does the default work out OK?)
* Analysis scripts can accept a parameter like `--output_prefix alternate_` to prefix all the output filenames with `alternate_`. That's handy to distinguish them from a previous run.
* You can run a single analysis script via

       python models/ecoli/analysis/cohort/growthDynamics.py

   or run one or more of them like this:

       runscripts/manual/analysisCohort.py growthDynamics transcriptFrequency
* The manual scripts don't compress any output files.
