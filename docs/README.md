# Whole Cell Model of E. coli

These are the docs for a variety of topics on the Whole Cell Model.

## setup


There are two alternative ways to set up to run the model.

Building a Docker Container Image takes one step -- one shell script.
The Container will be isolated from your computer's operating system, versions of Python, binary libraries, and everything else installed on your development computer.

Building a pyenv virtual environment takes a lot of installation steps that depend
on what's already installed. The result is isolated from other Python
virtual environments but the binary libraries might impact the rest of your computer.
On the other hand, this lets you edit, run, and debug without the complications of
Docker mechanics.

* **Docker setup**

   Install the [Docker Desktop software](https://www.docker.com/products/docker-desktop).

   **On Windows:** To use Docker on Windows, do it within
   [Windows Subsystem for Linux 2 (WSL2)](https://docs.microsoft.com/en-us/windows/wsl/about)
   which is a real Linux kernel that behaves the same as other Linux environments.
   Don't try to use VirtualBox.

   **NOTE:** Open Docker's Advanced Preferences and increase the memory allocation to 4GB.
   (The default allocation is 2GB which would make the model's Python code run out of
   memory, print "Killed", and stop with exit code 137.)

   Build and run the Docker Container Image locally like this:

   ```shell script
   cloud/build-containers-locally.sh
   docker run --name=wcm -it --rm wcm-code
   ```

   Or build a Docker Image using a Google Cloud Build server,
   pull the Image from the Google Cloud Package Registry,
   and run it locally
   (to do this you'll need to [set up the Google Cloud project](google-cloud.md)):

   ```shell script
   cloud/build.sh
   docker pull gcr.io/${PROJECT}/${USER}-wcm-code
   docker run --name=wcm -it --rm gcr.io/${PROJECT}/${USER}-wcm-code
   ```

   **NOTE:** It's better to build the Docker Image on Linux (e.g. in Google Cloud
   Build or WSL2) than on macOS where `build-containers-locally.sh` has to work around
   a Docker bug by disabling AVX2 instructions in OpenBLAS. That change somehow
   makes the computation produce different results and run a little slower.
   An Image built on Linux won't have the AVX2 workaround and yet it runs fine on macOS.
   See [xianyi/OpenBLAS#2244](https://github.com/xianyi/OpenBLAS/issues/2244#issuecomment-696510557).

   Once you have a Docker Image, you can run the model's Python programs inside the Container.
   (PyCharm Pro should support debugging into a Docker Container but we haven't tested that.)

   After changing the model's source code in the `wcEcoli/` directory, building an
   updated `wcm-code` Container Image via Google Cloud Build takes only a few minutes:

   ```shell script
   cloud/build-wcm.sh
   ```

   **TIP:** To preserve the model's output files after the Container exits,
   bind its output directory `/wcEcoli/out` to a host directory like `out/` by adding
   the option `-v $PWD/out:/wcEcoli/out`, where `$PWD` is the
   path to your cloned repo in the host computer.

   **NOTE:** `-v` needs absolute paths!

   You can share the entire `/wcEcoli` directory to also substitute the model's code
   inside the Container with the code in your host wcEcoli directory by changing
   that option to `-v $PWD:/wcEcoli`.

   **TIP:** When running in Docker on Linux or WSL2, the Container's output files
   will have the wrong host user and group ownership. You can change that by
   adding the `--user "$(id -u):$(id -g)"` option to the `docker run` command.
   That runs the process inside the Container as your host computer user and group so
   the files will be owned by you, but it adds other complications since it runs
   inside the Container without ownership of the existing files and directories.

* **pyenv setup**

  1. [Required development tools](dev-tools.md) to install the development tools including pyenv, gcc, make, and git, then

  1. [Creating the "pyenv" runtime environment](create-pyenv.md) to set up the Python runtime virtual environment for the model including binary libraries and Python packages.
`pyenv` virtual environments take more steps to build and depend on your computer's OS, but are lighter weight to run and easier for debugging.

   You can then run the model with this version of Python.

   If you have Anaconda installed, you might have to take Anaconda off the `$PATH` temporarily to run the Whole Cell Model.

   This approach takes a bunch of steps that vary depending on your operating system. It will run ≈25% faster than inside a Container and works with any Python debugger.

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
