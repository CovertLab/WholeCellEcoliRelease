# Surrogates

Light weight cell models that can be added to an environmental simulation.

# Uses

Surrogates are used for developing some feature of either the environment or the cell model that is focused on cell-environment interactions. The whole cell model can also be compared with these much simpler models to see what added predictive value it provides.

# Execution

Add surrogates into an environmental simulation using ```experiment``` from agent framework. This is the same way the whole cell model is added to an environmental simulation: 

    python -m environment.boot experiment --number N --type T
        
Here, ```T``` specifies the surrogate type, such as ```chemotaxis```, and ```N``` specifies the number of cells.

Timelines are an optional argument:

    python -m environment.boot experiment --number N --type T --timeline L

Here, ```L``` would be a string specifying events with time (seconds) and media_id: ```'0 minimal, 100 minimal_plus_amino_acids, 200 minimal'```

They can also be added to an already-running experiment with ```add```:

    python -m environment.boot add --type T

    
# Making New Surrogates

In order to make a new surrogate, create a new python file in this directory. It requires some interface with ```agent.inner```, including the functions: 
- ```run_incremental(run_until)```, for running the cell until the time specified by run_until. 
- ```generate_inner_update()```, for passing values such as environment_change to the environment.
- ```apply_outer_update(update)```, for passing updates from the environment to the cell.

A new surrogate also requires the addition of functions to ```environment.boot```:
- ```boot_SURROGATE(agent_id, agent_type, agent_config)```, for setting up the interface with ```agent.inner```
- ```initialize_SURROGATE(agent_id, agent_type, agent_config)```, for passing configuration data to the surrogate, calling the boot function, and starting the simulation.

Finally, it also requires an initializer:
- ```initializers['SURROGATE'] = initialize_SURROGATE```