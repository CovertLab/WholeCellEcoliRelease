# Development tools

This page explains how to install the needed development tools. You won't have to revisit this very often.


### First requirement: A package manager

**macOS:** You'll need the [Homebrew](https://brew.sh/) package manager to install software tools and libraries.
See the [Homebrew website](https://brew.sh/) about supported releases of macOS and how to install Homebrew.
There is a substitute for Homebrew on older releases of macOS but you might still run into SDK API
incompatibilities. We can't help you there short of suggesting that you update to the current
release of macOS.

**Linux:** You'll need to use `apt` or another package manager.

**Windows:** _We have not tested our software on Windows._ You could run our software within a Linux
virtual machine or go it alone on Windows. We won't be able to help you there.


## Required tools: gcc, make, git

### For macOS

You'll need Xcode's command line tools including the C compiler, make, and git.

Run this shell command to install _or update_ the command line tools:

```bash
xcode-select --install
```


### For Ubuntu

```bash
sudo apt install -y gcc make build-essential wget curl llvm
```

For Sherlock, the needed tools are already installed.

For other OSs, you'll need at least **gcc**, **git**, and **make**.



## Required tools: pyenv and virtualenv

_pyenv_ and _virtualenv_ are tools to install versions of Python and switch between virtual environments. A virtual environment contains an installation of Python and additional libraries. Separate projects should use separate virtual environments so they won't clash over these requirements.

1. Install `pyenv`, `pyenv-virtualenv`, `pyenv-virtualenvwrapper` using your local package manager, e.g. [homebrew](https://brew.sh/) on macOS. E.g.
   ```bash
   brew install pyenv pyenv-virtualenv pyenv-virtualenvwrapper
   ```

   See [pyenv Installation](https://github.com/pyenv/pyenv#installation) for various ways to install pyenv on various operating systems.

   Or try these commands on Ubuntu:

   ```bash
   git clone https://github.com/pyenv/pyenv.git ~/.pyenv
   echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bash_profile
   echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bash_profile
   echo -e 'if command -v pyenv 1>/dev/null 2>&1; then\n  eval "$(pyenv init -)"\nfi' >> ~/.bash_profile
   source ~/.bash_profile

   git clone https://github.com/pyenv/pyenv-virtualenv.git $(pyenv root)/plugins/pyenv-virtualenv
   echo 'eval "$(pyenv virtualenv-init -)"' >> ~/.bash_profile

   git clone https://github.com/pyenv/pyenv-virtualenvwrapper.git $(pyenv root)/plugins/pyenv-virtualenvwrapper
   source ~/.bash_profile
   ```

2. Set your shell login script (`~/.bash_profile` on Linux; `~/.profile` or `~/.bash_profile` on macOS, etc.) to initialize `pyenv` and optionally `pyenv-virtualenv` for each shell. To do this, follow the steps below or the more intricate instructions under "Add pyenv init to your shell" in [pyenv Installation](https://github.com/pyenv/pyenv#installation).
   - Example `~/.profile` or `~/.bash_profile` lines for macOS:

   ```bash
   export PYENV_ROOT=/usr/local/var/pyenv
   if which pyenv > /dev/null; then eval "$(pyenv init -)"; fi
   if which pyenv-virtualenv-init > /dev/null; then eval "$(pyenv virtualenv-init -)"; fi
   ## -- Do this *before* sourcing iterm2_shell_integration
   ```

   - Example `~/.bash_profile` lines for Stanford's "Sherlock" compute cluster, and note the `module load` command (more info at [setup_getting_started.md](https://github.com/CovertLab/ComputationalResources/blob/master/_sherlock/setup_getting_started.md), but that doc is out of date):

   ```bash
   module load wcEcoli/sherlock2

   export PYENV_ROOT="${PI_HOME}/pyenv"

   if [ -d "${PYENV_ROOT}" ]; then
       export PATH="${PYENV_ROOT}/bin:${PATH}"
       eval "$(pyenv init -)"
       eval "$(pyenv virtualenv-init -)"
   fi
   ```

3. Open a new shell so it runs the updated profile.

4. You'll need to put the project on the `PYTHONPATH` when working on it. Consider adding this to your profile _or_ creating a shell alias to do it when you work on wcEcoli:

   ```bash
   export PYTHONPATH="$HOME/wcEcoli:$PYTHONPATH"
   ```

5. **NOTE:** if you have a `~/.local/` directory, paths might not work properly with `pyenv` and you might receive error messages. ([TODO] What error messages? Delete the directory?)


## Additional requirements for macOS Mojave

**macOS Mojave (10.14+):**
  * **Update FUSE before upgrading to macOS Mojave** if you have FUSE.
  * After upgrading to Mojave, run `xcode-select --install` again.
  * On Mojave, follow the "you will also need to install the additional SDK headers" instructions on [https://github.com/pyenv/pyenv/wiki/Common-build-problems](https://github.com/pyenv/pyenv/wiki/Common-build-problems). In short:

   ```bash
   sudo installer -pkg /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg -target /
   ```


## Recommended tools

  * [PyCharm](https://www.jetbrains.com/pycharm/) or PyCharm Pro -- a very productive Integrated Development Environment (IDE)
  * [Sublime Text 3](https://www.sublimetext.com/) -- a slick code editor with many speed features
  * [GitHub Desktop app](https://desktop.github.com/) -- greases the skids for common git operations and lets you compose commit messages while reviewing the edits
  * [iTerm2](https://www.iterm2.com/) for macOS -- much more helpful than the stock Terminal app


### PyCharm setup

After [building the pyenv](docs/create-pyenv.md) and cloning the repo to a local directory, you can create a project in PyCharm.

* Select the project's Python interpreter: PyCharm > Preferences > Project: wcEcoli > Project Interpreter > Add > **Existing environment** > [navigate to something like `/usr/local/var/pyenv/versions/wcEcoli2/python`].
* Set Keyboard Shortcuts: Duplicate one of the available sets (e.g. "Mac OS X 10.5+"), then make changes to suit. A great change is to set Cmd-D (or Ctrl-D) for "Edit > Find > Add Selection for Next Occurrence". [This is like `find_under_expand` in Sublime Text. Cmd-G (or Ctrl-G) works like Sublime's `find_under_expand_skip`.]

[TODO] Notes on setting up code style, inspections, ...


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
