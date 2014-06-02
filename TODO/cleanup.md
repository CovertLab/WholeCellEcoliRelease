
Repository clean-up
-------------------

1. Top-level directories
    * Clean up TODO directory
    * fixtures - should these all go in out/fixtures?
    * runscripts - break into subdirs

2. wholecell/processes
    * Remove unneeded/outmoded processes
    * Clean up old code left in processes

3. wholecell/reconstruction
    * Move spreadsheet to a better location (kbEcoli repo?)
    * Move model_observations.md to a better location (top-level documentation directory?)
    * Move processes.txt to a better location (top-level documentation directory?)

4. wholecell/sim/hooks.py
    * Create a separate directory for hooks
    * Separate file per hook

5. wholecell/tests
    * Eliminate deprecation warnings (DONE?)
    * More/better tags

6. wholecell/utils
   * Eliminate alternative polymerize implementations
   * Remove configfile directory (need to fix Hudson scripts) (DONE)
   * Move some config options to constants (DONE)
   * Eliminate fixture manager (DONE)


Refactoring
-----------

1. wholecell/analysis
    * Functions for common operations (IO in particular)
    * Tests for analysis scripts that assert IO behavior
    * Separate directories for analysis scripts and analysis modules

2. wholecell/containers
    * Move saving/loading to bulk_objects_container from bulk_molecules
    * Rewrite Transcripts, Chromosome to use lists-of-arrays instead of single arrays
    * Rewrite UniqueObjectsContainer to better differentiate the global reference array, improved array expansion logic using numpy built-ins, and a global name parameter that allows for unique object IDs
    * Move the data loader into a new directory and rewrite to focus on the container instead of the state

3. wholecell/listeners
    * Saving/loading functionality

4. wholecell/loggers/shell.py
    * Implement many possible logged attributes and allow passing as a list of strings during simulation instantiation

5. wholecell/processes
    * Cache request calculations for use in evolveState (esp. elongation processes)

6. wholecell/reconstruction
    * Eliminate redundant/dependent entries from the KB object
    * Do something about validation
    * Represent arrays as sparse and use methods to access full arrays (created on demand)
    * Move more information out of KB into SQL
    * Move more information into KB from processes
    * Create unique molecule entries based on attribute lists + names of bulk molecules (as to inherit the proper masses)
    * Correct any misuse of unit'd values

7. wholecell/states
    * Standardize state/view interaction
    * Track mass in each state, by subtype/process/pre,post-evolveState/compartment

8. wholecell/utils
    * Do something about flex_t_fba_model.py (DONE - Derek is in charge)
    * Do something about linear_programming.py (DONE - see above)
    * Eliminate rand_stream and use numpy RNG objects (DONE)

9. Partitioning
    * Change bulk molecules partitioning to the new algorithm
    * Move partitioning calculations to the simulation level
    * Partitioner as a class
    * Coupling options (none, groups within/across states)
    * Priority options
    * Listener for partitioning (requests, satisfaction, etc.)
    * Better heuristics for degenerate cases (see code in UniqueMolecules)

10. Mass
    * Need to track more quantities, both for fitting and mass conservation
    * Unified update calculations for unique instances
