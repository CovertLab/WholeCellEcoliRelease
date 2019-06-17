# Whole Cell Model - *Escherichia coli*

This repository contains work to date on the [Covert Lab's](https://www.covert.stanford.edu/) Whole Cell Model for [*Escherichia coli*](https://en.wikipedia.org/wiki/Escherichia_coli), as well as some effort to create a framework for building whole cell models in general.

You can reach us at [AllenCenterCovertLab](mailto:allencentercovertlab@gmail.com).


## Setup

See [docs/README.md](docs/README.md) for docs on how to set up and run the model.

In short, there are two alternative ways to set up to run the model: inside a Docker container vs. in a carefully constructed `pyenv` virtual environment.


## Quick start

When running this code, prepare with these steps (the wcm-code Docker container already prepares this for you):

1. `cd` to the top level of your `wcEcoli` directory.
2. Set the `$PYTHONPATH`:

   ```bash
   export PYTHONPATH="$PWD:$PYTHONPATH"
   ```

3. In the `wcEcoli` directory, compile the Cython code:

   ```bash
   make clean compile
   ```


There are three ways to run the model:

   1. Use the manual runscripts.  
      They run each step directly, which is particularly handy to use with a
      debugger. But you're responsible for properly sequencing the steps.
   2. Queue up a Fireworks workflow, then run it.  
      This will automatically run all the steps, including parameter calculation,
      multiple cell simulations, and all the analysis plots.
   3. Use the multi-scale agent-based framework.  
      This runs multiple cells interactively on a simulated microscope slide.


## Using the manual runscripts

These scripts will run the parameter calculator (ParCa), cell simulations, and analysis steps directly, without any workflow. They're handy for development, e.g. running under the PyCharm debugger. But you're responsible for running the scripts in order and for re-running the ParCa after relevant code changes.

These scripts have command line interfaces built on `argparse`, so you can use shorter option names as long as they're unambiguous, and also one-letter forms so you can use `--cpus 8`, or `--cpu 8`, or `-c8`.

**NOTE:** _Use the `-h` or `--help` switch to get complete, up-to-date documentation on the command line options._ Below are just _some_ of the command line options.


To run the parameter calculator (ParCa), which is needed to prepare data for the simulation:
```bash
python runscripts/manual/runParca.py [-h] [--cpus CPUS] [--cached] [sim_outdir]
```

To simulate one or more cell generations with optional variants:

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

NOTE: The analysis plots currently look for multigen and variant output files
in the directories written by a Fireworks workflow, not the directories written
by `runSim.py`.

If you default the parameters, the analysis scripts will pick the latest simulation directory, the first variant, the first generation, and so on.
The list of PLOTs can include analysis class filenames like `aaCounts` (or `aaCounts.py`)
and analysis group TAGS like `CORE`. The default is to run the `CORE` group of plots that are
recommended for everyday development.

Set the environment variable `DEBUG_GC=1` if you want to check for Python memory
leaks when running the analysis plots.

The `--plot` (or `-p`) optional parameter lets you pick one or more specific plots to run. For example, to run two analysis plots on simulation variant #3 and put a filename prefix "v3_" on their output files (to distinguish them from other analysis runs):

```bash
python runscripts/manual/analysisCohort.py --plot compositionFitting.py figure2e.py --variant_index 3 --output_prefix v3_
```

There's another way run an individual analysis plot:

```bash
python models/ecoli/analysis/cohort/transcriptFrequency.py [-h] [--verbose] [-o OUTPUT_PREFIX] [-v VARIANT_INDEX] [sim_dir]
```

## Causality

After running a simulation, see the [CovertLab/causality](https://github.com/CovertLab/causality) repo for the Causality visualization tool to examine the model's causal links and simulation output correlations.


## Running an entire Fireworks workflow

See [wholecell/fireworks/README.md](wholecell/fireworks/README.md) for instructions to set up MongoDB as needed to run Fireworks.

The command line program `fw_queue.py` queues up a Fireworks workflow including parameter calculations, the simulation itself, and analysis plots.

The `fw_queue.py` source code begins with documentation on its _many options._
The options are set via environment variables. Below are a few usage examples.

But first, note that you can reset the Fireworks queue (if needed) via:

```bash
lpad reset
```

### Single simulation

To queue up a single simulation in Fireworks, including parameter calculations and analysis plots:

```bash
DESC="Example run of a single simulation." python runscripts/fireworks/fw_queue.py
```

The `DESC` text should be more descriptive than this so you can readily distinguish your runs.

### Multiple simulations

To queue multiple simulations, e.g. 4 simulations, each with a different initial seed:

```bash
DESC="Example run of multiple simulations." N_INIT_SIMS=4 python runscripts/fireworks/fw_queue.py
```

### Multiple generations

To queue multiple generations, e.g. 4 generations from a single mother cell:

```bash
DESC="Example run of multiple generations." N_GENS=4 python runscripts/fireworks/fw_queue.py
```

To queue multiple generations (in this case 3 generations) from multiple mother cells (in this case 2 mother cells:

```bash
DESC="Example run of multiple generations from multiple mother cells." N_GENS=3 N_INIT_SIMS=2 python runscripts/fireworks/fw_queue.py
```

### Shifting media conditions

To queue a simulation that switches between environments, use the "timeline" variant and give the range of indices (in this case from 1 to 1) specifying conditions defined in wcEcoli/environment/condition/timelines:

```bash
DESC="Example run of nutrient shifts." VARIANT="timeline" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 python runscripts/fireworks/fw_queue.py
```

### Using the cached sim data object

To use the cached sim data file, set the `CACHED_SIM_DATA` environment variable
(TODO: Explain what creates a cached sim data file):

```bash
DESC="Example run with cached sim data." CACHED_SIM_DATA=1 python runscripts/fireworks/fw_queue.py
```

### Using an interactive Sherlock node to run a Fireworks workflow

To run queued simulations on an interactive Sherlock node:

```bash
rlaunch rapidfire
```

You probably only want to do this if you're running or debugging a single simulation
(one initial seed, generation, and variant).

Don't do this on a Sherlock login node.

### Using the SLURM scheduler on Linux to run a Fireworks workflow

To run simulations on a Sherlock cluster (helpful when running more than one simulation):

```bash
qlaunch -r rapidfire --nlaunches infinite --sleep 5
```

The `qlaunch` command will run forever. Hit `Ctrl-C` to kill it once the console
logs shows that all the simulation and analysis steps have finished.

`qlaunch` is relatively lightweight, so it might work on a Sherlock login node.

`qlaunch` will create block directories with stdout and stderr from each Firework.  To troubleshoot any errors or just to see the output you would normally see from an interactive session, use the following commands to search the block directories for your desired fw_id:
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


## Multi-scale agent framework

See [How to run multi-scale agents](environment/README.md).
