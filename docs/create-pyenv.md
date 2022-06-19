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

The [requirements.txt](https://github.com/CovertLab/wcEcoli/blob/master/requirements.txt) file lists the Python packages needed for the Python runtime environment, with terse setup instructions.

This page goes through the Python environment setup steps in more detail and with more options, but `requirements.txt` gets updated more frequently since that file specifies the current package versions.

**NOTE**: While you can create virtual environments using
`virtualenv` or `venv` in place of `pyenv`, be sure to put the environment
*outside* the `wcEcoli/` directory. Otherwise `make clean` will break it!

**Sherlock:** Sherlock is the Stanford scientific computing cluster. Outside the Covert lab, just skip our Sherlock notes. Inside the lab, look in `$GROUP_HOME/downloads/` and `$GROUP_HOME/installation_notes/` for downloaded software packages and notes on recompiling them as needed to install new packages, new libraries, and new Python releases for the team.


## Prerequisites

* **Install** the software tools as described in [dev-tools](dev-tools.md), including
  * pyenv and pyenv-virtualenv
  * initializing pyenv in your shell profile
  * gcc or llvm
  * git
  * a programming editor such as PyCharm, Sublime Text, or Visual Studio Code
* **[Set up Git and GitHub](https://docs.github.com/en/get-started/quickstart/set-up-git)** including [Connecting to GitHub with SSH](https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh).
* **Clone the repo** [wcEcoli git](https://github.com/CovertLab/wcEcoli) to a local directory like `~/dev/wcEcoli/`. See [About remote repositories](https://docs.github.com/en/get-started/getting-started-with-git/about-remote-repositories). It's probably better to use the SSH URL `git@github.com:CovertLab/wcEcoli.git` but it's also possible to use an HTTPS URL.


## Install native libraries

1. Use your package manager to install the needed libraries
\[see the `requirements.txt` file for the latest list] or compile them from source.

   Most of this list comes from pyenv's requirements to install Python releases,
   so check [the pyenv wiki](https://github.com/pyenv/pyenv/wiki) for the latest
   list of libraries.

   We no longer recommend installing `openblas` library using a package manager or
   by compiling from source. Instead we let NumPy and SciPy install their
   own copy and let Aesara find NumPy's copy.

   You can optionally install `openblas` via package manager or source code,
   but be sure to get at least release v0.3.9.

   **On macOS**

   ```bash
   brew install glpk openssl readline swig suite-sparse xz
   ```

   **On Ubuntu**

   ```bash
   sudo apt install -y glpk-utils libglpk-dev glpk-doc libssl-dev libreadline-dev \
     libncurses5-dev libncursesw5-dev libffi-dev zlib1g-dev libbz2-dev xz-utils \
     libsqlite3-dev tk-dev
   ```

   For Ubuntu, you might also need to find and install the proprietary package `python-glpk`.


   Don't use apt-get to install `libopenblas-dev` until that package repository
   updates to a recent release of OpenBLAS such as the one that's embedded in
   numpy and scipy package wheels.
   On Ubuntu, use ```apt-cache policy libopenblas-dev``` to check the candidate version.
   To see which version of OpenBLAS is embedded in numpy, see
   [openblas_support.py](https://github.com/numpy/numpy/blob/main/tools/openblas_support.py)
   in the relevant numpy release tag.
   Recommendation: Let numpy and scipy install their embedded copies of OpenBLAS
   (see below).


   **On Sherlock**

   The needed packages are already installed. Set up your bash profile to locate the
   group environment modules, load the git and python modules, and initialize `pyenv`.
   You'll need these newer git modules since they use a compatible version of `libressl`.

   ```shell script
   export PI_HOME=$GROUP_HOME

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

   See `$GROUP_HOME/installation_notes/python3.txt`.

   If you need to update binary libraries like libressl, readline, or libffi,
   see their `$GROUP_HOME/installation_notes/*.txt` files.

   Each of these libraries and tools needs an _environment module_. We
   `module load` the module to make it accessible via environment variable paths
   like `CPPFLAGS`. See for example `$GROUP_HOME/modules/xz/5.2.5.lua`.

### On Ubuntu and other Linux environments

1. Install Python using `pyenv`, enabling Python to be loaded as a shared
   library.

   ```bash
   PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 3.8.7
   ```

### Anywhere else

Note: If running Parca gets an error message that "the loader can't load libpython"
or an error related to "shared objects" and "needing to compile with -fPIC", use
the above command that has `--enable-shared`.

1. Use `pyenv`.

   ```bash
   pyenv install 3.8.7 
   ```


## Create the "wcEcoli3" python virtual environment

**On Sherlock:** Skip this section and instead run `pyenv local wcEcoli3` to set up
your project directory unless you need to update or rebuild the team's shared
virtualenv.


1. Create a pyenv virtualenv and select it for your project directory.

   ```bash
   cd ~/dev/wcEcoli  # or wherever you cloned the `wcEcoli` project to
   pyenv virtualenv 3.8.7 wcEcoli3 && pyenv local wcEcoli3
   ```

1. Upgrade this virtual environment's installers.

   ```bash
   pip install --upgrade pip setuptools virtualenv virtualenvwrapper virtualenv-clone wheel
   ```

1. ***CONDITIONAL:*** Link numpy and scipy to a manually-installed OpenBLAS.

   **See Issue #931.** There are several degrees of freedom for installing
   OpenBLAS, numpy, and scipy that might change the computed results. We do not know
   how to set up environments to get consistent results across platforms.
   The simplest and fastest setup is to install numpy and scipy from binary "wheels"
   with their embedded copies of OpenBLAS using `pip install` **without** the
   `--no-binary numpy,scipy` option.
   Using the same BLAS implementation should help with portable reproducibility,
   and numpy's embedded copy is compiled with gcc/gfortran for every platform, so
   this is the recommended way to go.

   Still, there's a case for using a package manager to install OpenBLAS or
   compiling OpenBLAS from source code, then linking numpy and scipy to it as
   `cloud/docker/runtime/Dockerfile` _optionally_ does when building the
   wcm-runtime Docker Image.

   In that case, you need to locate the OpenBLAS library.
   * Brew on macOS installs OpenBLAS in `/usr/local/opt/openblas/`.
   * For other package managers, find out where they installed `lib/libopenblas*`.
   * Compiling OpenBLAS from source in Ubuntu goes into `/opt/OpenBLAS/` by default.
   * Compiling from source with
     `make FC=gfortran && make PREFIX=/XYZ install` installs it in that specified `/XYZ`
     PREFIX directory.
   * On Sherlock, it's installed in `$GROUP_HOME/downloads-sherlock2/compiled/openblas`.
     (Using an environment module to load the OpenBLAS when installing numpy and
     scipy works if the same environment module is loaded at runtime.)

   To link numpy and scipy to a manually-installed OpenBLAS, create a `~/.numpy-site.cfg` file pointing to
   it and remember to run `pip install <packages> --no-binary numpy,scipy` in the
   pip-install steps below.

   Copy these lines to `~/.numpy-site.cfg`, adjusting the paths as needed:

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
   pip install numpy==1.19.5  # see requirements.txt for the right version
   ```

   **NOTE:** If you installed OpenBLAS and created `~/.numpy-site.cfg`, use this command
   instead so pip will compile numpy from source code using `~/.numpy-site.cfg`:

   ```shell script
   pip install numpy==1.19.5 --no-binary numpy  # see requirements.txt for the right version
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

   The `LDFLAGS="-shared $LDFLAGS"` preamble fixes dozens of scipy build errors starting with  
   `In function _start (.text+0x20): undefined reference to main` and  
   `undefined reference to PyFloat_FromDouble`.

1. **Prerequisite for the following Python steps:** Set the `PYTHONPATH` environment variable, e.g.

   ```shell script
   export PYTHONPATH=$PWD
   ```

   or run the `ppath` alias recommended in [dev-tools.md](dev-tools.md).

1. Test the NumPy and SciPy installation.

      ```shell script
      python runscripts/debug/summarize_environment.py
      ```

      It should print entries like the ones below for numpy and scipy showing which
      OpenBLAS they're linked to. ```library_dirs = ['/usr/local/opt/openblas/lib']```
      reveals the path to source openblas while ```library_dirs = ['/usr/local/lib']```
      is shown for numpy's embedded openblas.

      ```
      lapack_opt_info:
          libraries = ['openblas', 'openblas']
          library_dirs = ['/usr/local/opt/openblas/lib']
          define_macros = [('HAVE_CBLAS', None)]
          language = c
      ```
      or this:

      ```
      lapack_opt_info:
         libraries = ['openblas', 'openblas']
         library_dirs = ['/usr/local/lib']
         language = c
         define_macros = [('HAVE_CBLAS', None)]
      ```

      (NumPy's embedded OpenBLAS library is in the virtualenv's
      `site-packages/numpy.libs/` on Linux, or
      `site-packages/numpy/.dylibs/` on macOS, or
      `site-packages/numpy/.libs/` on Windows.)


1. **Required:** Add the following line to your shell profile and run it in your current shell.
This gets more consistent results from OpenBLAS and it improves performance significantly,
especially when called from multiple processes.

    ```shell script
    export OPENBLAS_NUM_THREADS=1
    ```

1. Time the NumPy and SciPy libraries

    ```shell script
    runscripts/debug/time_libraries.sh
    ```

    (It might fail some timing expectations with a message like "AssertionError: 0.xxxxx not less than or equal to 0.4.")

1. Test Aesara:

    ```shell script
    python -c 'import aesara; print([aesara.config.blas.ldflags, aesara.config.device, aesara.config.floatX])'
    ```

   It should print something like this (with variations if you compiled OpenBLAS from source):

    `['-lblas', 'cpu', 'float64']`

1. Compile the project's native code.

   ```shell script
   make clean compile
   ```

   (Yes, it's expected to print deprecation warnings.)

1. Run the unit tests and runParca.py.

   ```shell script
   pytest

   ppath
   python runscripts/manual/runParca.py
   ```

   If the unit tests or runParca.py fail with an error message saying the loader can't load
   ...libpython... or an error related to "shared objects" and "needing to compile with -fPIC",
   that means you need to `--enable-shared` when installing python.
   Go back to that step, run

   ```shell script
   PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 3.8.7 --force
   ```

   then delete and recreate the virtualenv `wcEcoli3`.
   Delete it via the command `pyenv virtualenv-delete wcEcoli3` or `pyenv uninstall wcEcoli3`.

1. Remember to copy the git hooks from the repo (see [git hooks](../runscripts/git_hooks/README.md)) to your `.git` directory to maintain an up to date environment while doing development:

   ```
   cp runscripts/git_hooks/*[^.md] .git/hooks/
   ```

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

   ```shell script
   mkdir $SCRATCH/wcEcoli_out
   cd wcEcoli
   ln -s $SCRATCH/wcEcoli_out out
   ```

1. Create a symbolic link to a shared sim data cache directory in `$GROUP_SCRATCH` that should contain a copy of the newest sim data object (it should be updated by the daily build):

   ```shell script
   ln -s $GROUP_SCRATCH/wc_ecoli/cached cached
   ```
