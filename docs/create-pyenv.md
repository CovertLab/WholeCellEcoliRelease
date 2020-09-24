# Create the Python runtime environment

## NOTE: Docker is much easier

Before you begin all the steps below to install the Python runtime for the Whole Cell Model, consider how much easier it is to run it within a Docker container.

The script `cloud/build-containers-locally.sh` will build a Docker Image on
your computer if you have Docker Engine or Docker Desktop installed.
The script `cloud/build.sh` will do the same using Google Container Registry if
you have a Google Cloud project set up.

Either way, you can run the Whole Cell Model (WCM) in a Docker Container, isolated from your computer's operating system, Python versions, binary libraries, and everything else installed.

See [docs/README](README.md).


## Background

The [requirements.txt](https://github.com/CovertLab/wcEcoli/blob/master/requirements.txt) file contains the needed package list for the Python runtime environment and terse setup instructions.

This page goes through the Python environment setup steps in more detail and with more options.

**Prerequisites:** Install the software tools as described in [dev-tools](dev-tools.md). That page covers installing pyenv and pyenv-virtualenv, initializing them in your shell profile, installing a C compiler, and more.

**NOTE**: While it is possible to create a virtual environment with
`virtualenv` or `venv` in place of `pyenv`, be sure to put the environment
*outside* the `wcEcoli/` directory. Otherwise `make clean` will break it!

**Sherlock:** Sherlock is the Stanford scientific computing cluster. Outside the Covert lab, just skip our Sherlock notes. Inside the lab, look in `$PI_HOME/downloads/` and `$PI_HOME/installation_notes/` for downloaded software packages and notes on recompiling them as needed to install new packages, new libraries, and new Python releases for the team.

**See Issue #931.** There are several degrees of freedom for installing
OpenBLAS, numpy, and scipy which change the computed results. We do not know
how to set up environments to get consistent results across platforms.
The simplest and fastest setup is to install numpy and scipy from binary "wheels"
with their embedded copies of OpenBLAS. Still, there's a case for compiling
OpenBLAS from source code and linking numpy and scipy to it, as
`cloud/docker/runtime/Dockerfile` does for building the wcm-runtime Docker
Image.


## Install native libraries

1. Use your package manager to install the needed libraries
\[see the `requirements.txt` file for the latest list] or compile them from source.

   Theano will use the `openblas` library installed in this step.
   You can optionally install numpy and scipy to also use it.

   **On macOS**

   ```bash
   brew install glpk openssl readline swig suite-sparse xz openblas
   ```

   **On Ubuntu**

   ```bash
   sudo apt install -y glpk-utils libglpk-dev glpk-doc libssl-dev libreadline-dev \
     libncurses5-dev libncursesw5-dev libffi-dev zlib1g-dev libbz2-dev xz-utils \
     libsqlite3-dev tk-dev openblas
   ```

   For Ubuntu, you might also need to find and install the proprietary package `python-glpk`.

   Don't use apt-get to install `libopenblas-dev` until that package repository
   updates to a recent release like v0.3.9 (the version that's embedded in numpy
   and scipy).

   **On Sherlock**

   The needed packages are already installed. Set up your bash profile to locate the
   group environment modules, load the git and python modules, and initialize `pyenv`.
   You'll need these newer git modules since they use a compatible version of `libressl`.

   ```shell script
   ##### Add group-wide path settings #####
   if [ -f "${PI_HOME}/etc/bash_profile" ]; then
       . "${PI_HOME}/etc/bash_profile"
   fi

   module load git/2.27.0 git-lfs/2.11.
   module load wcEcoli/python3

   export PYENV_ROOT="${PI_HOME}/pyenv"

   if [ -d "${PYENV_ROOT}" ]; then
       export PATH="${PYENV_ROOT}/bin:${PATH}"
       eval "$(pyenv init -)"
       eval "$(pyenv virtualenv-init -)"
   fi
   ```

1. Optional: Download and install other packages according to their instructions or take a wait-and-see approach with them.

   * CPLEX from IBM (free for students)


## Install Python

### On Sherlock

1. Install Python 3 **in a shared pyenv for the team** if it needs updating.

   See `$PI_HOME/installation_notes/python3.txt`.

   If you need to update binary libraries like libressl, readline, or libffi,
   see their `$PI_HOME/installation_notes/*.txt` files.

   Each of these libraries and tools needs an _environment module_. We
   `module load` the module to make it accessible via environment variable paths
   like `CPPFLAGS`. See for example `$PI_HOME/modules/xz/5.2.5.lua`.


### Anywhere else

1. Install Python using `pyenv`.
pyenv lets you install and switch between multiple Python releases and multiple
"virtual environments", each with its own pip packages.

   ```bash
   PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 3.8.5
   ```


## Create the "wcEcoli3" python virtual environment

**On Sherlock:** Skip this section and instead run `pyenv local wcEcoli3` to set up
your project directory unless you need to update or rebuild the team's shared
virtualenv.


1. Create a pyenv virtualenv and select it for your project directory.

   ```bash
   cd ~/dev/wcEcoli  # or wherever you cloned the `wcEcoli` project to
   pyenv virtualenv 3.8.5 wcEcoli3 && pyenv local wcEcoli3
   ```

1. Upgrade this virtual environment's installers.

   ```bash
   pip install --upgrade pip setuptools virtualenv virtualenvwrapper virtualenv-clone wheel
   ```

1. ***CONDITIONAL:*** Link numpy and scipy to a manually-installed OpenBLAS.

   **See Issue #931.** There are several degrees of freedom for installing
   OpenBLAS, numpy, and scipy which change the computed results. We do not know
   how to set up environments to get consistent results across platforms.
   The simplest and fastest setup is to install numpy and scipy from binary _wheels_
   with their embedded copies of OpenBLAS using `pip install` **without** the
   `--no-binary numpy,scipy` option. Still, there's a case for compiling
   OpenBLAS from source code and linking numpy and scipy to it, as
   `cloud/docker/runtime/Dockerfile` does for building the wcm-runtime Docker
   Image.

   Where is OpenBLAS installed?
   * Brew on macOS installs OpenBLAS in `/usr/local/opt/openblas/`.
   * For other package managers, find out where they install it.
   * When compiling OpenBLAS from source,
     `make FC=gfortran && make PREFIX=/XYZ install` installs it in that `/XYZ`
     PREFIX directory, which defaults to `/opt/OpenBLAS`.
   * On Sherlock, it's installed in `$PI_HOME/downloads-sherlock2/compiled/openblas`.
     (Using an environment module for OpenBLAS only works if it's loaded at runtime.)

   To link numpy and scipy to a manually-installed OpenBLAS, create a `~/.numpy-site.cfg` file pointing to
   it (and remember to run `pip install <packages> --no-binary numpy,scipy` in the
   pip-install steps below), e.g.:

      ```
      [openblas]
      libraries = openblas
      library_dirs = /usr/local/opt/openblas/lib
      include_dirs = /usr/local/opt/openblas/include
      runtime_library_dirs = /usr/local/opt/openblas/lib
      ```

1. Install NumPy.

   Install numpy before installing `scipy` and `stochastic-arrow` to avoid
   installation errors.

   ```shell script
   pip install numpy==1.19.2  # see requirements.txt for the right version
   ```

   **NOTE:** If you installed OpenBLAS and created `~/.numpy-site.cfg`, use this command
   instead so pip will compile numpy from source code using `~/.numpy-site.cfg`:

   ```shell script
   pip install numpy==1.19.2 --no-binary numpy  # see requirements.txt for the right version
   ```

1. Install the packages listed in `requirements.txt`.

   ```shell script
   pip install -r requirements.txt && pyenv rehash
   ```

   **NOTE:** If you installed OpenBLAS and created `~/.numpy-site.cfg`, use this command
   instead:

   ```shell script
   LDFLAGS="-shared $LDFLAGS" pip install -r requirements.txt --no-binary numpy,scipy && pyenv rehash
   ```

   The `LDFLAGS="-shared $LDFLAGS"` preamble fixes dozens of scipy build
   errors starting with  
   `In function _start (.text+0x20): undefined reference to main` and  
   `undefined reference to PyFloat_FromDouble`.

1. Test the NumPy and SciPy installation.

      ```bash
      python runscripts/debug/summarize_environment.py
      ```

      It should print entries like this for numpy and scipy showing which
      OpenBLAS they're linked to:

      ```
      lapack_opt_info:
          libraries = ['openblas', 'openblas']
          library_dirs = ['/usr/local/opt/openblas/lib']
          define_macros = [('HAVE_CBLAS', None)]
          language = c
      ```

1. **(Now required)** Add the following line to your bash profile and run it in your current shell.
This gets more consistent results from OpenBLAS and it improves performance significantly,
especially when called from multiple processes.

    ```
    export OPENBLAS_NUM_THREADS=1
    ```

1. Time the NumPy and SciPy libraries

    ```bash
    runscripts/debug/time_libraries.sh
    ```

1. Test Theano:

      ```bash
      python -c 'import theano; print(theano.config.blas.ldflags)'
      ```

   It should print something like

      `-lblas`

   or

      `-L/usr/local/opt/openblas/lib -lopenblas -lopenblas`

1. Compile the project's native code.

   ```bash
   make clean compile
   ```

   (Yes, it's expected to print deprecation warnings.)

1. Run the unit tests.

   ```bash
   export PYTHONPATH=$PWD
   pytest
   ```

   If the unit tests fail with an error message saying the loader can't load
   ...libpython..., that means you need to `--enable-shared` when installing python.
   Go back to that step, run

   ```bash
   PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 3.8.5 --force
   ```

   then delete and recreate the virtualenv `wcEcoli3`.

1. If you're using PyCharm, be sure to select the project's Python interpreter so PyCharm understands the version
of Python and its installed libraries. This enables code completion, usage documentation
in context, visual debugging, warnings about code problems, click-through to library
source code, etc.

   > PyCharm >  
   > Preferences >  
   > Project: wcEcoli >  
   > Project Interpreter >  
   > gear ⚙️ >  
   > Add... >  
   > Virtualenv Environment >  
   > Existing environment >  
   > Interpreter >  
   > [run `pyenv which python` in a shell to find the python location, something
   > like `/usr/local/var/pyenv/versions/wcEcoli3/python`, and paste that path
   > into the text box or navigate there].

## Sherlock SCRATCH directory setup

1. Make sure the model's output goes to the `$SCRATCH` filesystem (which is larger) rather than SHERLOCK HOME.

   ```bash
   mkdir $SCRATCH/wcEcoli_out
   cd wcEcoli
   ln -s $SCRATCH/wcEcoli_out out
   ```

1. Create a symbolic link to a shared sim data cache directory on `$PI_SCRATCH` that should contain a copy of the newest sim data object (it should be updated by the daily build):

   ```bash
   ln -s $PI_SCRATCH/wc_ecoli/cached cached
   ```
