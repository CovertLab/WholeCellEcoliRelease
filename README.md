# Whole Cell Model - *Escherichia coli*

**Notice:** This repository contains a **release snapshot** of the [Covert Lab's](https://www.covert.stanford.edu/) Whole Cell Model for [*Escherichia coli*](https://en.wikipedia.org/wiki/Escherichia_coli). In contrast, our working repository is under continuous development so please contact us before embarking on any changes that you want to contribute. We do not plan to merge Pull Requests into this repository except documentation and installation fixes.

You can reach us at [AllenCenterCovertLab](mailto:allencentercovertlab@gmail.com).


## Setup

See [docs/README.md](docs/README.md) for docs on how to set up and run the model. As a shortcut, the core setup instructions are in:

  1. [Required development tools](docs/dev-tools.md) to install the basic tools, then
  2. [Creating the "pyenv" runtime environment](docs/create-pyenv.md) to set up the Python runtime environment for the model.


## Quick start

Before running this code, remember to:

1. `cd` to the top level of your cloned `wcEcoli` directory, and
2. set the `$PYTHONPATH`

   ```bash
   export PYTHONPATH="/path/to/wcEcoli:$PYTHONPATH"
   ```

In the `wcEcoli` directory, compile the Cython code:

```bash
make clean compile
```


There are two ways to run the model:

   1. Use the manual runscripts.
   2. Queue up a Fireworks workflow, then run it.


## Using the manual runscripts

These scripts will run the parameter calculator (misnamed "fitter"; it does not do machine learning), simulation, and analysis steps directly, without any workflow. They're handy for development, e.g. running under the PyCharm debugger. But when you use them, you're responsible for running the scripts in order and for re-running the parameter calculator after relevant code changes.

These scripts have command line interfaces built on Python's `argparse`, so you can use argparse features like shorter option names as long as they're unambiguous. Many options also have one-letter forms like `-c8` as short for `--cpus 8`.

**NOTE:** _Use the `-h` or `--help` switch to get complete, up-to-date documentation on the command line options._ Below are just _some_ of the command line options.


To run all the parameter calculation steps:
```bash
python runscripts/manual/runFitter.py [-h] [--cpus CPUS] [sim_outdir]
```

To do a simple simulation run:

```bash
python runscripts/manual/runSim.py [-h] [--variant VARIANT_TYPE FIRST_INDEX LAST_INDEX] [--generations GENERATIONS] [--seed SEED] [sim_dir]
```

To run all the analysis plots on the simulation output in a given `sim_dir`:

```bash
python runscripts/manual/analysisCohort.py [-h] [--plot PLOT [PLOT ...]] [--cpus CPUS] [--output_prefix OUTPUT_PREFIX] [--variant_index VARIANT_INDEX] [sim_dir]

python runscripts/manual/analysisMultigen.py [-h] [--plot PLOT [PLOT ...]] [--cpus CPUS] [--output_prefix OUTPUT_PREFIX] [--variant_index VARIANT_INDEX] [--seed SEED] [sim_dir]

python runscripts/manual/analysisSingle.py [-h] [--plot PLOT [PLOT ...]] [--cpus CPUS] [--output_prefix OUTPUT_PREFIX] [--variant_index VARIANT_INDEX] [--seed SEED] [--generation GENERATION] [--daughter DAUGHTER] [sim_dir]

python runscripts/manual/analysisVariant.py [-h] [--plot PLOT [PLOT ...]] [--cpus CPUS] [--output_prefix OUTPUT_PREFIX] [sim_dir]
```

If you default the parameters, the script will pick the latest simulation directory, the first variant, the first generation simulation, and so on.

Set the environment variable `DEBUG_GC=1` to check for Python memory leaks in the analysis scripts.

The `--plot` (or `-p`) optional parameter lets you pick one or more specific plots to run. For example, to run two analysis scripts on simulation variant #3 and put a filename prefix "v3_" on their output files (to distinguish them from other analysis runs):

```bash
python runscripts/manual/analysisCohort.py --plot compositionFitting.py figure2e.py --variant_index 3 --output_prefix v3_
```

You can also run an individual analysis script directly:

```bash
python models/ecoli/analysis/cohort/transcriptFrequency.py [-h] [--verbose] [-o OUTPUT_PREFIX] [-v VARIANT_INDEX] [sim_dir]
```


## Running an entire Fireworks workflow

See [wholecell/fireworks/README.md](wholecell/fireworks/README.md) for instructions to set up MongoDB et al to run Fireworks.

The command line program `fw_queue.py` queues up a Fireworks workflow including parameter calculations, the simulation itself, and lots of analysis plots.

The `fw_queue.py` source code begins with documentation on its many options. Below are a few usage examples.

But first, note that you can reset the Fireworks queue (if needed) via:

```bash
lpad reset
```

### Single simulation

To queue up a single simulation in Fireworks, including parameter calculations and analysis plots:

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

### Using an interactive session to run a Fireworks workflow

After queuing your simulation(s), to run them in an interactive session, run:

```bash
rlaunch rapidfire
```

You probably only want to do this if you're running/debugging a single simulation.

### Using the SLURM workload scheduler on Linux to run a Fireworks workflow

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
