(1) Biomass objective function 


(2) Protein decay and the number of ribosomes required to double the cell protein fraction
The number of ribosomes needed to double the protein fraction of the cell based on expression, protein degradation, and dilution was calculated in the fitter based on the following equation:
```
nRibosomesNeeded = np.sum(
		monomerLengths / kb.ribosomeElongationRate * (
			np.log(2) / kb.cellCycleLen + kb.monomerData["degRate"]
			) * monomersView.counts()
		).to('dimensionless').magnitude
```
