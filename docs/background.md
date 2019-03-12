# Background on the Whole Cell Model


## Source code locations

* The knowledge base input information: `wcEcoli/reconstruction/ecoli/flat`.
* The code that uses it to calculate parameters: `wcEcoli/reconstruction/ecoli/`.
* The code that loads it into an object: `knowledge_base_raw.py`.
* Simulation processes: `wcEcoli/models/ecoli/processes/`.
* The cell's initial conditions: `wcEcoli/models/ecoli/sim/initial\_conditions.py`.


## Processes

Each _process_ models one part of the cell’s function, e.g. RNA polymerase elongation. They are modeled separately (modular), run in short time steps, and the results from each time step are integrated between processes before initiating the next time step. Other processes include translation, elongation, transciption, and metabolism.

Each process has three entry points:
* _initialize_: called only once at the beginning of a simulation. Get needed parameters from the knowledge base, get views of bulk and unique molecules (bulk molecules are “indistinguishable” from each other, e.g. inactive RNAP molecules, unique molecules can be distinguished from each other, e.g. active RNAP molecules are each assigned to a location on the genome), create a view so that you can get counts, change counts, and change properties.
* _calculateRequest_: called at the beginning of each timestep. Request the resources that you want for that timestep (don’t request all unless you are certain that another process doesn’t need this resource as well, don’t forget about metabolism).
* _evolveState_: called after resources are allocated at each timestep. Perform the process, update counts, and update masses (mass must be conserved between steps).


## Listeners

If you want to save new values to analyze after the simulation, you must write it out via a _listener_. Listeners log data during the simulation. To add a listener to the code, add the file to `wcEcoli/models/ecoli/listeners/` and add the listener to `_listenerClasses` in `wcEcoli/models/ecoli/sim/simulation.py`.


## States, bulk molecules, and unique molecules

[TODO] Introduce States, bulk molecules, and unique molecules.


## Design factors

* One implicit modeling design goal is that no phenomena be modeled in more than one place over a given time interval -- exactly one place, if the model is complete.

* We have to be careful about degrees of representation. E.g. genes with known differential expression but no associated transcription factor do not currently change their expression during an environmental shift.

* A reproducible/testable way to show what is or isn't represented is via validation against some expected behavior or withheld data set.
