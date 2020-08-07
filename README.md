# Whole Cell Model - *Escherichia coli*

**Notice:** This repository contains a **release snapshot** of the [Covert Lab's](https://www.covert.stanford.edu/) Whole Cell Model for [*Escherichia coli*](https://en.wikipedia.org/wiki/Escherichia_coli). In contrast, our working repository is under continuous development so please contact us before embarking on any changes that you want to contribute. We do **not** plan to merge Pull Requests into this repository except documentation and installation fixes.

You can reach us at [AllenCenterCovertLab](mailto:allencentercovertlab@gmail.com).


## Setup

See [docs/README.md](docs/README.md) for more info on setting up and running the model.

In short, there are two alternative setups to run the model: inside a Docker container vs. in a manually constructed `pyenv` virtual environment.


## Quick start

When running this code, prepare with these steps (the wcm-code Docker container already prepares this for you):

1. `cd` to the top level of your `wcEcoli` directory.
2. Set the `$PYTHONPATH`:

   ```bash
   export PYTHONPATH="$PWD"
   ```

3. In the `wcEcoli` directory, compile the Cython code:

   ```bash
   make clean compile
   ```


There are two ways to run the model:

   1. Use the manual runscripts.
   2. Queue up a FireWorks workflow, then run it.


## Using the manual runscripts

These scripts will run the parameter calculator (misnamed "fitter"; it does not do machine learning), simulation, and analysis steps directly, without any workflow. They're handy for development, e.g. running under the PyCharm debugger. But you're responsible for running the steps in order and for re-running the relevant steps after changes.
Or you can delete the output directory and run the steps from the beginning.

These scripts have command line interfaces built on Python's `argparse`, so you can use argparse features like shorter option names as long as they're unambiguous. Many options also have one-letter forms like `-c8` as short for `--cpus 8`.

**NOTE:** _Use the `-h` or `--help` switch to get complete, up-to-date documentation on the command line options._ Below are just _some_ of the command line options.


To run the parameter calculator (ParCa), which is needed to prepare input data for the simulation:
```bash
python runscripts/manual/runFitter.py [-h] [--cpus CPUS] [sim_outdir]
```

To simulate one or more cell generations with optional variants:

```bash
python runscripts/manual/runSim.py [-h] [--variant VARIANT_TYPE FIRST_INDEX LAST_INDEX] [--generations GENERATIONS] [--seed SEED] [sim_dir]
```

To run all the analysis plots on the simulation output in a given `sim_dir`:

```bash
python runscripts/manual/analysisVariant.py [-h] [--plot PLOT [PLOT ...]] [--cpus CPUS] [sim_dir]

python runscripts/manual/analysisCohort.py [-h] [--plot PLOT [PLOT ...]] [--cpus CPUS] [--variant_index VARIANT_INDEX] [sim_dir]

python runscripts/manual/analysisMultigen.py [-h] [--plot PLOT [PLOT ...]] [--cpus CPUS] [--variant_index VARIANT_INDEX] [--seed SEED] [sim_dir]

python runscripts/manual/analysisSingle.py [-h] [--plot PLOT [PLOT ...]] [--cpus CPUS] [--variant_index VARIANT_INDEX] [--seed SEED] [--generation GENERATION] [--daughter DAUGHTER] [sim_dir]
```

If you default the parameters, an analysis script will pick the latest simulation directory, the first variant, the first generation simulation, and so on.

To get the full set of plot outputs after running multiple variants, seeds, and/or
generations, run:
* `analysisVariant` once;
* `analysisCohort` on each variant;
* `analysisMultigen` on each combination of variant × seed;
* `analysisSingle` on each combination of variant × seed × generation.

Set the environment variable `DEBUG_GC=1` to check for Python memory leaks in the analysis scripts.

The `--plot` (or `-p`) optional parameter lets you run one or more specific scripts from a category. For example, to run two analysis scripts on simulation variant #3 and put a filename prefix "v3_" on their output files (to distinguish them from other analysis runs):

```bash
python runscripts/manual/analysisCohort.py --plot compositionFitting.py figure2e.py --variant_index 3 --output_prefix v3_
```


## Running an entire FireWorks workflow

See [wholecell/fireworks/README.md](wholecell/fireworks/README.md) for instructions to set up MongoDB to run FireWorks.

The command line program `fw_queue.py` queues up a FireWorks workflow (using the
MongoDB server configuration info in `my_launchpad.yaml`) including parameter
calculations, the simulation itself, and the full set of analysis plots.

The `fw_queue.py` source code begins with documentation on its many options. Below are a few usage examples.

But first, note that you can reset the FireWorks queue (if needed) via:

```bash
lpad reset
```

### Single simulation

To queue up a single simulation in FireWorks, including parameter calculations and analysis plots:

```bash
DESC="Example run of a single simulation." python runscripts/fw_queue.py
```

The `DESC` text should be descriptive so you can readily distinguish your runs from each other.

### Multiple simulations

To queue multiple simulations (in this case 4 simulations):

```bash
DESC="Example run of multiple simulations." N_INIT_SIMS=4 python runscripts/fw_queue.py
```

### Multiple generations

To queue multiple generations (in this case 4 generations) from a single mother cell:

```bash
DESC="Example run of multiple generations." N_GENS=4 python runscripts/fw_queue.py
```

To queue multiple generations from multiple mother cells (in this case 3 generations each from 2 mother cells):

```bash
DESC="Example run of multiple generations from multiple mother cells." N_GENS=3 N_INIT_SIMS=2 python runscripts/fw_queue.py
```

### Shifting nutrient conditions

To queue a simulation that switches between environments, use the "nutrient_time_series" variant, and give the range of indices (in this case from 1 to 1) specifying conditions defined in wcEcoli/reconstruction/ecoli/flat/condition/timeseries:

```bash
DESC="Example run of nutrient shifts." VARIANT="nutrient_time_series" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 python runscripts/fw_queue.py
```

### Using an interactive session to run a FireWorks workflow

After queuing your simulation(s), to run them in an interactive session, run:

```bash
rlaunch rapidfire
```

You probably only want to do this if you're running/debugging a single simulation.

### Using the SLURM workload scheduler on Linux to run a FireWorks workflow

To run simulations using the cluster (you'll probably want to do this if you're running more than one simulation and/or more than one generation), run:

```bash
qlaunch -r rapidfire --nlaunches infinite --sleep 5
```

This command will run forever until you `Ctrl-C` to kill it once you see that all the output and analysis files have been generated.

`qlaunch` will create block directories with stdout and stderr from each firework.  To troubleshoot any errors or just to see the output you would normally see from an interactive session, use the following commands to search the block directories for your desired fw_id:
```bash
./runscripts/fw_info.sh out 1
./runscripts/fw_info.sh error 1
```
This will display the stdout and stderr from the execution of a firework with fw_id of 1.


### Output

The output is stored as a time-stamped sub-directory of the `out` directory, for example `out/20180703.215222.029168__multi-sim/`, where `DESC="multi-sim"` was one of the arguments to `fw_queue.py`.

Within this directory, there is a `metadata` sub-directory which stores the git revision information as well as the description provided by the `DESC` variable, a `kb` sub-directory which stores kb objects (after the simulations and analysis are done the objects are compressed using bzip2), and sub-directories (maybe only a single sub-directory) containing different variants (e.g., gene knockouts or other perturbations).

Within variant sub-directories, there are `N_INIT_SIMS` (which defaults to 1) numbered sub-directories such as `000000` corresponding to "family trees".

Within each "family tree" sub-directory are generation sub-directories such as `generation_000000`.

Within each generation sub-directory are numbered individual simulation directories that contain `simOut` (for simulation data) and `plotOut` (for plots) sub-directories.

A family tree for 3 generations showing the relationship between numbered individual simulations is shown here:

```
gen 0:                 0
gen 1:         0               1
gen 2:     0       1       2       3
```
