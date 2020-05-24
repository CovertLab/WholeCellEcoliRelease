# reflect

Reflect on python objects.

## object_tree

Examine and transform the structure of a python object into a dictionary.

### usage

First, import `object_tree`:

```python
import runscripts.reflect.object_tree as o
```

Then, get the object you wish to examine. In our case we are going to take a look at `sim_data`:

```python
import cPickle
sim_data = cPickle.load(open('out/manual/kb/simData.cPickle', "rb"))
```

Now that we have our object, we can transform it into a nested dictionary:

```python
sim_tree = o.object_tree(sim_data)
```

You can pretty-print the result, optionally to an output stream:

```python
import pprint as pp
pp.pprint(sim_tree)

with open('user/dump.txt', 'w') as f:
    pp.pprint(sim_tree, f)
```

The `sim_tree` dictionary that results has the same structure as `sim_data`, but all objects that aren't leaf types are now dictionaries with a `!type` key containing the original type. The full list of default leaf types are defined in the `leaf_types` tuple in `object_tree.py`.

If you want, you can get some information about the structure as it is being explored. The second argument is `path`, which establishes the name of the root path, and the third argument is `debug`, which takes two options. If supplied with `'ALL'`, it will print out every path encountered during the traversal. If debug is set to `'CALLABLE'`, it will print out only the callable leaves.

```python
sim_tree = o.object_tree(sim_data, 'sim_data', 'CALLABLE')
```

besides building the `sim_tree` dictionary, will print out all instance variables that contain callables:

```
sim_data.process.two_component_system.derivatives: <function derivatives at 0x11f9771b8>
sim_data.process.two_component_system.derivativesJacobian: <function derivativesJacobian at 0x11f977230>
sim_data.process.two_component_system.derivativesParca: <function derivatives at 0x11f977488>
sim_data.process.two_component_system.derivativesParcaJacobian: <function derivativesJacobian at 0x11f977500>
sim_data.process.equilibrium.derivatives: <function derivatives at 0x11f977398>
sim_data.process.equilibrium.derivativesJacobian: <function derivativesJacobian at 0x11f977578>
```

## model_inspection

Inspect implementations in the model based on raw_data and sim_data and save to tsv files.

### usage

Specify inputs (raw_data and sim_data) and output directory:

```bash
runscripts/reflect/model_inspection.py -r path/to/raw_data.cPickle -s path/to/sim_data.cPickle -o output_dir/
```

If inputs are not specified, the script will run the parca to produce the output.  If output directory is not specified, the files will be saved in the script directory.

Check options with:

```bash
runscripts/reflect/model_inspection.py -h
```
