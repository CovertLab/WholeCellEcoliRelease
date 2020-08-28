# Running the model

See the top level [README](../README.md) for general instructions and [docs/README](README.md) for details on setup and running
and tradeoffs between different ways to run the model.


## FireWorks on Google Cloud

See [How to run the Whole Cell Model on the Google Cloud Platform](google-cloud.md)
for instructions to run a FireWorks workflow of cell simulations and analysis plots in Google Cloud or your local computer.

## Fireworks on Sherlock

See [Setting up to run FireWorks](wholecell/fireworks/README.md) for instructions to run a FireWorks workflow of cell simulations and analysis plots on Sherlock or your local computer. Start with `fw_queue.py`.

**NOTE:** If you get this error message connecting to an older MongoDB server
(such as on `mlab.com`, not the one in Google Cloud) with the pymongo library 3.9.0+
that wcEcoli uses in Python 3:

> pymongo.errors.OperationFailure: This MongoDB deployment does not support retryable writes. Please add retryWrites=false to your connection string.

the fix is to add these lines to your `my_launchpad.yaml` file:

```
mongoclient_kwargs:
  retryWrites: false
```

## Using the Manual Runscripts

[TODO] Rewrite this. Meanwhile, see the top level [README](../README.md).

You can accomplish largely the same thing that running via Fireworks does by running:

1. `runParca.py` -- creates the `kb/*` sim data files for the control case
2. `runSim.py` -- creates the modified modified sim data files for requested variants and runs the sims
3. `analysisParca.py`, `analysisCohort.py`, `analysisMultigen.py`, `analysisSingle.py`, `analysisVariant.py` -- generates the output plots; to generate all the plots, you must run some of these multiple times

When using the manual runscripts:

* It runs directly without FireWorks or MongoDB, but it doesn't distribute multiple sims and analysis plots to run in parallel across multiple computers.
* You can run it in a debugger.
* You can run just the parts you want and rerun them as needed. The manual scripts don't automate dependency management of data files flowing between steps such as simulation outputs as analysis plot inputs. It's on you to rerun code if things change, runSim before analysis. (Some analysis scripts get confused if the sim runs are more varied than expected. See [https://github.com/CovertLab/wcEcoli/issues/199](#199).)
* The scripts have command line help and are more rigorous at parsing args than `fw_queue.py` is at reading environment variables. Argparse provides helpful features like the ability to use any unambiguous option name prefix so you can use `--var` in place of `--variant_index`.
* The scripts do smart defaulting. runParca defaults to creating the sim dir `out/manual/` rather than a timestamped directory name. The others default to finding the most likely sim dir, i.e. the latest timestamped subdirectory of `out/` or else the alphabetically first subdirectory of `out/`. You can pass in a directory name that's absolute, or relative to `wcEcoli/`, or relative to `wcEcoli/out/`.
* Analysis scripts can accept a parameter like `--variant_index 1` to look for a subdirectory like `condition_000001/`.
* You can run a single analysis script via

       python models/ecoli/analysis/cohort/growthDynamics.py

   or run one or more of them like this:

       runscripts/manual/analysisCohort.py growthDynamics transcriptFrequency
* The manual scripts don't compress the output files.
