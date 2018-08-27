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
sim_data = cPickle.load(open('out/manual/kb/simData_Fit_1.cPickle', "rb"))
```

Now that we have our object, we can transform it into a nested dictionary:

```python
sim_tree = o.object_tree(sim_data)
```

The `sim_tree` dictionary that results has the same structure as `sim_data`, but all objects that aren't leaf types are now dictionaries with a `!type` key containing the original type. The full list of default leaf types are defined in the `leaf_types` tuple in `object_tree.py`.

If you want, you can get some information about the structure as it is being explored. The second argument is `path`, which establishes the name of the root path, and the third argument is `debug`, which takes two options. If supplied with `'ALL'`, it will print out every path encountered during the traversal. If debug is set to `'CALLABLE'`, it will print out only those leaves which hold functions or methods.

```python
sim_tree = o.object_tree(sim_data, 'sim_data', 'CALLABLE')
```

besides building the `sim_tree` dictionary, will print out all leaf methods:

```python
sim_data._addConditionData
sim_data._addHardCodedAttributes
sim_data.constants._buildConstants
sim_data.external_state._buildEnvironment
sim_data.external_state._getNutrientData
sim_data.getter._buildAllMasses
sim_data.getter._buildLocations
sim_data.getter.getLocation
sim_data.getter.getMass
sim_data.growthRateParameters.getDnaCriticalMass
sim_data.growthRateParameters.getFractionActiveRibosome
sim_data.growthRateParameters.getFractionActiveRnap
sim_data.growthRateParameters.getFractionIncreaseRnapProteins
sim_data.growthRateParameters.getRibosomeElongationRate
sim_data.growthRateParameters.getRnapElongationRate
sim_data.initialize
sim_data.internal_state._buildBulkMolecules
sim_data.internal_state._buildCompartments
sim_data.internal_state._buildUniqueMolecules
sim_data.internal_state.bulkMolecules.addToBulkState
sim_data.internal_state.uniqueMolecules.addToUniqueState
sim_data.mass._buildCDPeriod
sim_data.mass._buildConstants
sim_data.mass._buildDependentConstants
sim_data.mass._buildSubMasses
sim_data.mass._buildTrnaData
sim_data.mass._calculateGrowthRateDependentDnaMass
sim_data.mass._chromosomeSequenceMass
sim_data.mass._clipTau_d
sim_data.mass._getFitParameters
sim_data.mass._getTrnaAbundanceAtGrowthRate
sim_data.mass.getAvgCellDryMass
sim_data.mass.getBiomassAsConcentrations
sim_data.mass.getFractionMass
sim_data.mass.getMassFraction
sim_data.mass.getTrnaDistribution
sim_data.moleculeGroups._buildMoleculeGroups
sim_data.moleculeIds._buildMoleculeIds
sim_data.process.complexation._buildStoichMatrixMonomers
sim_data.process.complexation._findColumn
sim_data.process.complexation._findRow
sim_data.process.complexation._moleculeRecursiveSearch
sim_data.process.complexation.getMonomers
sim_data.process.complexation.massBalance
sim_data.process.complexation.massMatrix
sim_data.process.complexation.stoichMatrix
sim_data.process.complexation.stoichMatrixMonomers
sim_data.process.equilibrium._findColumn
sim_data.process.equilibrium._findRow
sim_data.process.equilibrium._makeDerivative
sim_data.process.equilibrium._makeMatrices
sim_data.process.equilibrium._moleculeRecursiveSearch
sim_data.process.equilibrium._populateDerivativeAndJacobian
sim_data.process.equilibrium._solveSS
sim_data.process.equilibrium.derivatives
sim_data.process.equilibrium.derivativesJacobian
sim_data.process.equilibrium.fluxesAndMoleculesToSS
sim_data.process.equilibrium.getFwdRate
sim_data.process.equilibrium.getMetabolite
sim_data.process.equilibrium.getMetaboliteCoeff
sim_data.process.equilibrium.getMonomers
sim_data.process.equilibrium.getRevRate
sim_data.process.equilibrium.getUnbound
sim_data.process.equilibrium.massBalance
sim_data.process.equilibrium.massMatrix
sim_data.process.equilibrium.setFwdRate
sim_data.process.equilibrium.setRevRate
sim_data.process.equilibrium.stoichMatrix
sim_data.process.equilibrium.stoichMatrixMonomers
sim_data.process.metabolism._buildBiomass
sim_data.process.metabolism._buildMetabolism
sim_data.process.metabolism.concentrationUpdates._addMoleculeAmounts
sim_data.process.metabolism.concentrationUpdates._isNutrientExchangePresent
sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients
sim_data.process.metabolism.exchangeConstraints
sim_data.process.metabolism.getKineticConstraints
sim_data.process.replication._buildGeneData
sim_data.process.replication._buildReplication
sim_data.process.replication._buildSequence
sim_data.process.replication._reverseComplement
sim_data.process.rna_decay._buildRnaDecayData
sim_data.process.rna_decay.kmLossFunction
sim_data.process.transcription._buildRnaData
sim_data.process.transcription._buildTranscription
sim_data.process.transcription_regulation._buildLookups
sim_data.process.transcription_regulation.pPromoterBoundSKd
sim_data.process.transcription_regulation.pPromoterBoundTF
sim_data.process.translation._buildMonomerData
sim_data.process.translation._buildTranslation
sim_data.process.translation._buildTranslationEfficiency
sim_data.process.two_component_system._buildComplexToMonomer
sim_data.process.two_component_system._makeDerivative
sim_data.process.two_component_system._makeDerivativeFitter
sim_data.process.two_component_system._populateDerivativeAndJacobian
sim_data.process.two_component_system.derivatives
sim_data.process.two_component_system.derivativesFitter
sim_data.process.two_component_system.derivativesFitterJacobian
sim_data.process.two_component_system.derivativesJacobian
sim_data.process.two_component_system.getMonomers
sim_data.process.two_component_system.getReactionName
sim_data.process.two_component_system.makeDependencyMatrix
sim_data.process.two_component_system.massBalance
sim_data.process.two_component_system.massMatrix
sim_data.process.two_component_system.moleculesToNextTimeStep
sim_data.process.two_component_system.moleculesToSS
sim_data.process.two_component_system.stoichMatrix
sim_data.process.two_component_system.stoichMatrixMonomers
sim_data.relation._buildMonomerIndexToRnaMapping
sim_data.relation._buildRnaIndexToMonomerMapping
```