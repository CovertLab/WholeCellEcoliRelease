
Repository clean-up
-------------------

1. Top-level directories
    * Clean up TODO directory (DONE)
    * fixtures - eliminate? move them to out/fixtures?
    * runscripts - break into subdirs (DONE)

2. wholecell/processes (DONE)
    * Remove unneeded/outmoded processes (DONE)
    * Clean up old code left in processes (DONE)

3. wholecell/reconstruction (DONE)
    * Move spreadsheet to a better location (kbEcoli repo?) (DONE - moved to docs directory for now)
    * Move model_observations.md to a better location (top-level documentation directory?) (DONE)
    * Move processes.txt to a better location (top-level documentation directory?) (DONE)

4. wholecell/sim/hooks.py (DONE)
    * Create a separate directory for hooks (DONE)
    * Separate file per hook (DONE)

5. wholecell/tests (DONE)
    * Eliminate deprecation warnings (DONE?)
    * More/better tags (SKIPPING)

6. wholecell/utils (DONE)
   * Eliminate alternative polymerize implementations (DONE)
   * Remove configfile directory (need to fix Hudson scripts) (DONE)
   * Move some config options to constants (DONE)
   * Eliminate fixture manager (DONE)

7. Repository (DONE)
   * Remove/resolve floating branches (DONE)


Refactoring
-----------

1. wholecell/analysis
    * Functions for common operations (IO in particular)
    * Tests for analysis scripts that assert IO behavior
    * Separate directories for analysis scripts and analysis modules

2. wholecell/containers
    * Move saving/loading to bulk_objects_container from bulk_molecules
    * Save values related to requests/allocations/usages in an appropriate listener (PARTIALLY COMPLETE)
    * Rewrite Transcripts, Chromosome to use lists-of-arrays instead of single arrays (SKIPPING)
    * Rewrite UniqueObjectsContainer to better differentiate the global reference array (DONE), improved array expansion logic using numpy built-ins (DONE), and a global name parameter that allows for unique object IDs (DONE)
    * Move the data loader into a new directory and rewrite to focus on the container instead of the state

3. wholecell/listeners
    * Saving/loading functionality

4. wholecell/loggers/shell.py (DONE)
    * Implement many possible logged attributes and allow passing as a list of strings during simulation instantiation (have listeners register logged value formatting?) (DONE)

5. wholecell/processes (DONE)
    * Cache request calculations for use in evolveState (esp. elongation processes) (DONE) (but not used)

6. wholecell/reconstruction
    * Eliminate redundant/dependent entries from the KB object
    * Do something about validation (DONE - removed)
    * Represent arrays as sparse and use methods to access full arrays (created on demand) (DONE)
    * Move more information out of KB into SQL
    * Move more information into KB from processes
    * Create unique molecule entries based on attribute lists + names of bulk molecules (as to inherit the proper masses) (partially done)
    * Correct any misuse of unit'd values

7. wholecell/states
    * Standardize state/view interaction
    * Track mass in each state, by subtype/process/pre,post-evolveState/compartment (DONE except for compartments)

8. wholecell/utils (DONE)
    * Do something about flex_t_fba_model.py (DONE - Derek is in charge)
    * Do something about linear_programming.py (DONE - see above)
    * Eliminate rand_stream and use numpy RNG objects (DONE)

9. Partitioning
    * Change bulk molecules partitioning to the new algorithm
    * Move partitioning calculations to the simulation level
    * Partitioner as a class
    * Coupling options (none, groups within/across states)
    * Priority options (currently implemented for bulk molecules only)
    * Listener for partitioning (requests, satisfaction, etc.)
    * Better heuristics for degenerate cases (see code in UniqueMolecules)

10. Mass
    * Need to track more quantities, both for fitting and mass conservation (DONE)
    * Unified update calculations for unique instances

11. Simulation instantiation
    * Better ways to call limited sets of objects (i.e. states, processes, listeners)

Fitter clean up
---------------
1. Build functions for repetitive fitter code