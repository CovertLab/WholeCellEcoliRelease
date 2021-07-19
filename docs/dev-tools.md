# Development tools

This page explains how to install the needed development tools. You won't have to revisit this very often.


## First requirement: A package manager

**macOS:** You'll need the [Homebrew](https://brew.sh/) package manager to install software tools and libraries.
See the [Homebrew website](https://brew.sh/) about supported releases of macOS and how to install Homebrew.

**Linux:** You'll need `apt` or `snap` or another package manager.

**Windows:** _We have not tested the Whole Cell Model on Windows._ You could try Windows Subsystem for Linux (WSL 2),
or Docker on WSL 2, or a Linux virtual machine, or blaze the trail on Windows.


## Required tools: gcc, make, git

### On macOS

You'll need Xcode's command line tools including the C compiler, make, and git.

Run this shell command to install _or update_ the command line tools:

```shell script
xcode-select --install
```


### On Ubuntu

```shell script
sudo apt install -y gcc make build-essential wget curl llvm
```

### On Sherlock

(Here we're assuming you already have a login account for the Sherlock scientific computing cluster. If not,
you might not need to bother with Sherlock except to view Jenkins CI build logs.)

The needed tools are already installed for the group.
Look in `$PI_HOME/downloads/`, `$PI_HOME/installation_notes/`, and `$PI_HOME/modules/`.

* Just add to your shell login script as described below.


### On other OSs

You'll need at least `gcc`, `git`, and `make`.
Probably `cmake` and `llvm` as well.





## Required tools: pyenv and pyenv-virtualenv

`pyenv` and `virtualenv` are tools to install versions of Python and switch
between virtual environments, each with its own selection of Python and libraries.
Separate projects need separate virtual environments to avoid clashes over
required libraries and library versions.

1. Install `pyenv`, `pyenv-virtualenv`, `pyenv-virtualenvwrapper` using your local
   package manager, e.g. [homebrew](https://brew.sh/) on macOS.
   **\[On Sherlock, pyenv is already installed. Skip to the next step.]** 

   Do see [pyenv Installation](https://github.com/pyenv/pyenv#installation) for
   up-to-date details to install pyenv on various operating systems.

   Example commands on MacOS:

   ```shell script
   brew install pyenv pyenv-virtualenv pyenv-virtualenvwrapper
   ```

   Example commands on Linux:

   ```shell script
   git clone https://github.com/pyenv/pyenv.git ~/.pyenv
   git clone https://github.com/pyenv/pyenv-virtualenv.git $(pyenv root)/plugins/pyenv-virtualenv
   git clone https://github.com/pyenv/pyenv-virtualenvwrapper.git $(pyenv root)/plugins/pyenv-virtualenvwrapper
   ```

1. Configure your shell's environment (`~/.profile` or `~/.bash_profile` for bash,
   or `~/.zshrc` for zsh, etc.) to initialize `pyenv` and optionally `pyenv-virtualenv`
   per the latest instructions in
   [pyenv Installation](https://github.com/pyenv/pyenv#basic-github-checkout)
   under **Configure your shell's environment for Pyenv** _and_ under
   **Add pyenv into your shell**.

   - Example `~/.zshrc` lines for macOS:

     ```shell script
     if which pyenv > /dev/null; then
         eval "$(pyenv init --path)"
         eval "$(pyenv init -)"
     fi
     if which pyenv-virtualenv-init > /dev/null; then
         eval "$(pyenv virtualenv-init -)"
     fi
     ```

   - Example bash `~/.profile` or `~/.bash_profile` lines for macOS:

     ```shell script
     export PYENV_ROOT=/usr/local/var/pyenv
     if which pyenv > /dev/null; then
         eval "$(pyenv init --path)"
         eval "$(pyenv init -)"
     fi
     if which pyenv-virtualenv-init > /dev/null; then
         eval "$(pyenv virtualenv-init -)"
     fi
     ## ^^^ Do this before sourcing iterm2_shell_integration
     ```

   - Example bash `~/.profile` or `~/.bash_profile` lines for Linux:

     ```shell script
     export PYENV_ROOT="$HOME/.pyenv"
     export PATH="$PYENV_ROOT/bin:$PATH"
     if which pyenv > /dev/null; then
         eval "$(pyenv init --path)"
         eval "$(pyenv init -)"
     fi
     if which pyenv-virtualenv-init > /dev/null; then
         eval "$(pyenv virtualenv-init -)"
     fi
     ```

   - Example `~/.bash_profile` lines for Sherlock:

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

1. Open a new shell so it runs the updated profile.

   - On macOS, If you use normally use an account that does not have
     write access to `/usr/local/` (this is sometimes done to protect
     your Homebrew installation) you may see this error:

     ```
     mkdir: /usr/local/var/pyenv/shims: Permission denied
     mkdir: /usr/local/var/pyenv/versions: Permission denied
     ```

     In that case, change the `PYENV_ROOT` setting in your shell profile and
     open a new shell:

     ```shell script
     export PYENV_ROOT="$HOME/.pyenv"
     ```

1. You'll need to put the project source code root on the `PYTHONPATH` when working on it. Consider adding this to your profile:

   ```shell script
   export PYTHONPATH="$HOME/wcEcoli"
   ```

   _or_ create a shell alias **and run it when you work on wcEcoli or any other Python project**:

   ```shell script
   alias ppath='export PYTHONPATH=$PWD'
   ```

1. **NOTE:** if you have a `~/.local/` directory, paths might not work properly with `pyenv` and you might receive error messages. ([TODO] What error messages? Delete the directory?)


### Additional requirements for macOS Mojave and later

**macOS Mojave (10.14+):**
  * **If you have FUSE installed, update it before upgrading to macOS Mojave.**
  * After upgrading to Mojave, run `xcode-select --install` again.
  * On Mojave, follow the "you will also need to install the additional SDK headers" instructions on [https://github.com/pyenv/pyenv/wiki/Common-build-problems](https://github.com/pyenv/pyenv/wiki/Common-build-problems). In short:

   ```shell script
   sudo installer -pkg /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg -target /
   ```


## Recommended tools

  * [PyCharm](https://www.jetbrains.com/pycharm/) or PyCharm Pro -- a very productive Integrated Development Environment (IDE). PyCharm Professional Edition is free for students and academic staff members.
  * [Sublime Text](https://www.sublimetext.com/) -- a slick code editor with many speed features
  * [Visual Studio Code](https://code.visualstudio.com/) -- a slick code editor with IDE plug-ins
  * [GitHub Desktop app](https://desktop.github.com/) -- greases the skids for common git operations and lets you compose commit messages while reviewing the edits
  * [iTerm2](https://www.iterm2.com/) for macOS -- much more helpful than the stock Terminal app


### Shell tips

Define shell aliases and environment variable settings in a file such as `$HOME/bin/shell-aliases.sh` and "source" it in your shell profile file. This way, you can edit `shell-aliases.sh` and `source` it without having to start a new shell.

Example `.profile.sh`:

```shell script
export PATH=$PATH:$HOME/bin:/usr/local/sbin

[[ -f $HOME/bin/shell-aliases.sh ]] && source $HOME/bin/shell-aliases.sh

export PYENV_ROOT=/usr/local/var/pyenv
if which pyenv > /dev/null; then eval "$(pyenv init -)"; fi

# If you're using macOS iTerm2.
# See https://iterm2.com/documentation-shell-integration.html
test -e "${HOME}/.iterm2_shell_integration.bash" && source "${HOME}/.iterm2_shell_integration.bash"
```

Example `bin/shell-aliases.sh`:

```shell script
export OPENBLAS_NUM_THREADS=1
alias ppath='export PYTHONPATH=$PWD'

# If you have a symlink from `bin/subl@` -> `/Applications/Sublime Text.app/Contents/SharedSupport/bin/subl`
export EDITOR='subl -w'

export LESS="--ignore-case -R"
export MORE="--ignore-case --quit-if-one-screen"
alias la='ls -AF'
alias ll='ls -lhF'
alias ls='ls -F'
alias l.='ls -d .*'
alias lt='ls -lhFt'

alias cdw='cd ~/dev/wcEcoli'
alias cdwo='cd ~/dev/wcEcoli/out'

alias df="df -h"
alias du="du -h"
alias grep='grep --color'

# Append to the history file on exit so multiple tabs get their history saved.
shopt -s histappend
```


### PyCharm setup

After [building the pyenv](docs/create-pyenv.md) and cloning the repo to a local directory, you can create a project in PyCharm.
wcEcoli has a project in source control.

* **Be sure to select the project's Python interpreter so PyCharm understands the version
of Python and its installed libraries.** This enables code completion, usage documentation
in context, visual debugging, warnings about code problems, click-through to library
source code, and many other features for working with Python code.

  PyCharm >  
  Preferences >  
  Project: wcEcoli >  
  Project Interpreter >  
  gear ⚙️ >  
  Add... >  
  Virtualenv Environment >  
  Existing environment >  
  Interpreter >  
  [run `pyenv which python` in a shell to find the python location, something
  like `/usr/local/var/pyenv/versions/wcEcoli3/python`, and paste that path into
  the text box or navigate there].

* Set Keyboard Shortcuts: Duplicate a keymap (e.g. "Mac OS X 10.5+"), then make changes to suit.

  Recommended:

  * `Cmd-D` (or Ctrl-D) "Edit > Find > Add Selection for Next Occurrence".  
    This is like `find_under_expand` in Sublime Text. Cmd-G (or Ctrl-G) works like Sublime's `find_under_expand_skip`.
  * `Cmd-[` and `Cmd-]`: Navigate Back and Forward, like in a web browser.
  * `F2`: Next Bookmark, `Shift-F2`: Previous Bookmark, `Cmd-F2`: Toggle Bookmarks.

See [this PyCharm debugger tutorial video](https://youtu.be/QJtWxm12Eo0) or
[this PyCharm debugger tutorial webpage](https://www.jetbrains.com/help/pycharm/debugging-your-first-python-application.html).

#### Mismatched package names

PyCharm's inspector gives false-positive error messages like:

> "Package containing module 'Bio' is not listed in project requirements"

because it assumes imported module names (`Bio`) match their PyPI package names
(`biopython`) except for a specific list of exceptions. The built-in exceptions
list doesn't include `stochastic-arrow`, `biopython`, etc., but we can add them.

**The fix:** Each time you install a PyCharm release, run the following shell command,
then restart PyCharm:  

> \[NOTE: So far, this script only works on macOS and maybe only for PyCharm Pro
installed by JetBrains Toolbox. It won't work with snap-installed PyCharm since
snap installs on a readonly file system. To support other installations,
enhance the script for where `python.jar` or `python-ce.jar` gets installed.]

    runscripts/tools/augment-pycharm-package-list.sh

See [You Track issue #PY-27985](https://youtrack.jetbrains.com/issue/PY-27985).


### Great PyCharm features to know

* Install the "Markdown support" plug-in to help with editing `*.md` files.
* Cmd-click (or Ctrl-click) a name to jump to its definition.
* Navigate > Class, Navigate > File, Navigate > Symbol -- jump to any class, file, or symbol defined in the project or its libraries. This supports "fuzzy matching" where you can type characters in the name to narrow down the list.
* Edit > Find > Find Usages -- find all references to the selected symbol.
* Double-press Shift -- search everywhere in the project.
* Refactor -- commands to do small code refactorings like rename a function or change its calling signature.
* Hit `F1` to get documentation on an identifier's definition.


### iTerm2 Tips

**iTerm2** is a macOS app that [greatly improves on the stock terminal app](https://www.iterm2.com/features.html).

Tips (also see [the FAQ](https://www.iterm2.com/faq.html)):

* If you configure it to save & load preferences to a Dropbox folder, you don't have to do much when switching to a new Mac.
* tmux Integration lets you make and adjust window split panes much more easily than typing tmux commands.
* [Shell Integration](https://www.iterm2.com/documentation-shell-integration.html) is very handy, but the regular setup didn't work quite right on Sherlock with the pyenv virtualenv shell prompt. So for Sherlock, just set up the "Triggers" as documented on that page. The "Prompt Detected" trigger is probably the most useful part since it lets you jump between shell prompts in the terminal output.

Example "Default" profile configuration for "Keys":
* Option up/down arrows: scroll one line up/down
* Shift left/right arrows: send hex code `0x02` or `0x06` (respectively) to move the cursor
* Control left/right arrows: Send `^[b` or `^[f` (respectively; `^[` is ESC) to move the cursor by words (assuming your `.inputrc` file uses EMACS style readline editing, which is the default)
* Option left/right arrows: Same as Control left/right arrows
* Command left/right arrows: send hex code `0x01` or `0x05` (respectively) to move the cursor to the start/end of line
* Left Option (alt) key: Esc+
* Right Option (alt) key: Normal

Example "Default" profile configuration for "Triggers":
* `@login.sherlock.stanford.edu's password:` --> Open Password Manager
* `Enter passphrase for key '/home/users/\w+/.ssh/id_rsa':` --> Open Password Manager
* `^\[(\w+)@(sh-[\w-]+) login! ([^\]]+)]\$ ` --> Prompt Detected
* ^ You can use the same regex for Report User & Host `\1@\2` and for Report Directory `\3`

Example overall configuration for "Keys":
* Control Tab, Control-Shift Tab: Next Tab, Previous Tab
* Shift up/down arrows: Scroll One Line Up/Down
* Command up/down arrow: Scroll to Top/End
* Control up/down arrows, Page Up/Down, Command Page Up/Down, Shift Page Up/Down: Scroll One Page Up/Down
* Command Home/End: Scroll to Top/End
