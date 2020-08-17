# Whole Cell Model of E. coli

These are the docs for a variety of topics on the Whole Cell Model.

## setup


There are two alternative ways to set up to run the model:

* **Docker setup**

   Install the [Docker Desktop software](https://www.docker.com/products/docker-desktop) then pull our Docker container image from the Google Cloud Package Registry `gcr.io`, or build it locally using the `cloud/build-containers-locally.sh` or `cloud/build.sh` shell scripts.

   A Docker container image builds automatically and is isolated from your computer's operating system, versions of Python, binary libraries, and everything else installed on your development computer.

   You can then run the model inside the Docker container.
   (PyCharm Pro should support debugging into a Docker container but we haven't tested that.)

   You can share a local directory with the code inside the Docker container, either just the model's output directory `/wcEcoli/out` in order to preserve its output files outside the container:

   ```docker run --name=wcm -v <LOCAL_WCECOLI>/out:/wcEcoli/out --user "$(id -u):$(id -g)" -it --rm wcm-code```

   (where `<LOCAL_WCECOLI>` denotes the local path to your cloned repo),
   or share the entire `/wcEcoli` directory to also substitute the model's code inside the container with the code in your local wcEcoli directory.

   ```docker run --name=wcm -v <LOCAL_WCECOLI>:/wcEcoli --user "$(id -u):$(id -g)" -it --rm wcm-code```

   The `--user "$(id -u):$(id -g)"` option runs the model inside the container as your user and group in the host computer so the output files will be owned by you.

   **NOTE:** If you encounter memory issues while using Docker Desktop (the default allocated memory is 2GB) causing the simulation to get killed midway, click the Docker icon > Preferences > Advanced > adjust memory to 4GB.

   **NOTE:** Docker Desktop for Windows is not currently compatible with VirtualBox.  If you use VirtualBox, try installing the legacy [Docker Toolbox](https://github.com/docker/toolbox/releases) instead.  You may also need to adjust the memory allocated to the VirtualBox VM (named 'default') that gets created.  In VirtualBox, select the 'default' VM and under system, change the base memory from 1 GB to 4 GB. 

* **pyenv setup**

  1. [Required development tools](dev-tools.md) to install the development tools including pyenv, gcc, make, and git, then

  1. [Creating the "pyenv" runtime environment](create-pyenv.md) to set up the Python runtime virtual environment for the model including binary libraries and Python packages.
`pyenv` virtual environments take more steps to build and depend on your computer's OS, but are lighter weight to run and easier for debugging.

   You can then run the model with this version of Python.

   If you have Anaconda installed, you might have to take Anaconda off the `$PATH` temporarily to run the Whole Cell Model.

   This approach takes a bunch of steps that vary depending on your operating system. It will run ≈25% faster than inside a container and works with any Python debugger.

  See:

   * [Required development tools](dev-tools.md) -- installation and tips
   * [Creating the "pyenv" runtime environment](create-pyenv.md)
   * [Setting up to run FireWorks](../wholecell/fireworks/README.md) -- needed only to run a FireWorks workflow of cell simulations and analysis plots
   * [Set up Zookeeper and Kafka](../agent/README.md) -- needed only for multi-scale agents

* **Also**

  After setting up the environment, copy the git hooks from the repo (see [git hooks](../runscripts/git_hooks/README.md)) to your `.git` directory to maintain an up to date environment while doing development:

  ```
  cp runscripts/git_hooks/*[^.md] .git/hooks/
  ```

## running

* [How to run the Whole Cell Model](run.md) (actually, the top level [README](../README.md) is more informative)
* [How to run the Whole Cell Model on the Google Cloud Platform](google-cloud.md)
* [How to run the Causality visualization tool](https://github.com/CovertLab/causality)
* [How to run multi-scale agents](../environment/README.md)

## development

* [Background on the model](background.md)
* [Coding style guide](style-guide.md)

## relevant papers

* [Simultaneous cross-evaluation of heterogeneous _E. coli_ datasets via mechanistic simulation](https://science.sciencemag.org/content/369/6502/eaav3751.full), _Science_, 24 July 2020
* [A Whole-Cell Computational Model Predicts Phenotype from Genotype](https://www.cell.com/cell/abstract/S0092-8674(12)00776-3), _Cell_, July 20, 2012

## dissertations
* _Computational Simulations of Whole Cells: Strategies for Framework Design and Model Parameterization_, John Mason
* _Development and Application of Whole-Cell Computational Models for Science and Engineering_, Jonathan Ross Karr
* _Toward a Whole-Cell Model of Escherichia coli_, Derek Macklin
* _Towards a Whole-Cell Model of Growth Rate and Cell Size Control in Escherichia coli_, Nicholas Ruggero
* _Transcriptional Regulation in Escherichia coli: A Systems Biology Approach_, Markus Covert
