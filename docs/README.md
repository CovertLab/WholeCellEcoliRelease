# Whole Cell Model of E. coli

These are the docs for a variety of topics on the Whole Cell Model.

## setup

There are two alternative ways to set up to run the model:

1. **Docker Container:** Install the
   [Docker Desktop software](https://www.docker.com/products/docker-desktop)
   and run the Docker Image as a Container.

   Pull the full docker image from the package registry in this repo:
   ```shell script
   docker pull docker.pkg.github.com/covertlab/wholecellecolirelease/wcm-full:latest
   ```

   You can then run the wcEcoli model inside the Container like this:

   ```shell script
   docker run --name=wcm -it --rm docker.pkg.github.com/covertlab/wholecellecolirelease/wcm-full
   ```

   The `-it` option starts an interactive shell.
   Alternatively, you can supply a shell command to run.

   Another useful option is `--rm`, which asks Docker to remove the Container on
   exit so you don't have to remember to delete old Containers.

   You can mount your local directory `wcEcoli/out/` into the Container to preserve the
   program's output files when the Container exits - just be sure to provide a full path
   to `out/` (eg. `$PWD/out`), not just a relative path from your current directory:

   ```shell script
   docker run --name=wcm -v $PWD/out:/wcEcoli/out -it docker.pkg.github.com/covertlab/wholecellecolirelease/wcm-full
   ```

   In this case, the output files will be owned by root. You can fix
   that by adding the option `--user "$(id -u):$(id -g)"` to run the process
   inside the Container as your user and group from the host computer so the
   output files will be owned by you, but that adds complications. E.g.
   the process inside the Container won't have a user profile and won't own the
   `wcEcoli/` directory.

   **NOTE:** If you encounter memory issues while using Docker Desktop (the default allocated memory is 2GB) and the simulation processes get killed midway, click the Docker icon > Preferences > Advanced > adjust memory to 4GB.

   **NOTE:** When setting up Docker Desktop for Windows, it is best to [use the WSL2 backend](https://docs.docker.com/docker-for-windows/wsl/). If not, you might run into compatability issues when trying to run Docker, especially if you have VirtualBox installed, which could require using the legacy [Docker Toolbox](https://github.com/docker/toolbox/releases).

   Inside the Container you can then run commands like these:

   ```shell script
   python runscripts/manual/runFitter.py
   python runscripts/manual/runSim.py
   python runscripts/manual/analysisSingle.py
   ```

   **Tip:** Eventually, you'll want to delete the Docker Image. Refer to the
   commands `docker image prune`, `docker image ls`, and `docker image rm`.

   **Tip:** You can build your own Docker image instead of the one provided using these steps:

   ```shell script
   cd $YOUR_CODE_PROJECTS_DIR/wcEcoli  # or wherever you cloned the wcEcoli project to
   cloud/build-containers.sh
   ```

2. **Python virtual environment:** Follow [Required development tools](dev-tools.md) to install the development tools including pyenv, gcc, make, and git, then follow [Creating the "pyenv" runtime environment](create-pyenv.md) to set up the Python runtime virtual environment for the model including binary libraries and Python packages.

   You can then run wcEcoli in this Python virtualenv.

   This approach takes many careful steps that vary depending on your operating
   system but it will run noticeably faster than inside a Docker Container.
   The native libraries and compilers will not be isolated from the rest of your
   computer but the virtualenv will be isolated from other Python environments.
   However if you have Anaconda installed, you might have to
   take it off the `$PATH` temporarily to run Python in this virtualenv.

   * [Required development tools](dev-tools.md) -- installation and tips
   * [Creating the "pyenv" runtime environment](create-pyenv.md)
   * [Setting up to run FireWorks](../wholecell/fireworks/README.md) -- needed only to run a FireWorks workflow of cell simulations and analysis plots

## running

* [How to run the Whole Cell Model](run.md); actually the top level [README](../README.md) is more informative

## development

* [Background on the model](background.md)
* [Coding style guide](style-guide.md)

## relevant papers

* [Simultaneous cross-evaluation of heterogeneous E. coli datasets via mechanistic simulation](https://science.sciencemag.org/content/369/6502/eaav3751.full), _Science_, 24 July 2020
* [A Whole-Cell Computational Model Predicts Phenotype from Genotype](https://www.cell.com/cell/abstract/S0092-8674(12)00776-3), _Cell_, July 20, 2012

## dissertations
* _Computational Simulations of Whole Cells: Strategies for Framework Design and Model Parameterization_, John Mason
* _Development and Application of Whole-Cell Computational Models for Science and Engineering_, Jonathan Ross Karr
* _Toward a Whole-Cell Model of Escherichia coli_, Derek Macklin
* _Towards a Whole-Cell Model of Growth Rate and Cell Size Control in Escherichia coli_, Nicholas Ruggero
* _Transcriptional Regulation in Escherichia coli: A Systems Biology Approach_, Markus Covert
