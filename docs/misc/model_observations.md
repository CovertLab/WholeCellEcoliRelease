(1) <b>Biomass objective function</b>
Biomass compositoin how does that vary from cell to cell. Using a fixed biomass composition for a single cell leads to issues for us.

(2) <b>How does protein expression at the beginning of a cell cycle correlate with protein expression at the end of the cycle?</b>
A protein highly expressed at the beginning should still be highly expressed at the end of the cell cycle (?).


(3) <b>Protein decay and the number of ribosomes required to double the cell protein fraction</b><br>
TODO: UPDATE BECAUSE THIS IS NOW FIXED by adding half lives to protein counts fitting<br>
The number of ribosomes needed to double the protein fraction of the cell based on expression, protein degradation, and dilution was calculated in the fitter based on the following equation:

```
nRibosomesNeeded = np.sum(
		monomerLengths / kb.ribosomeElongationRate * (
			np.log(2) / kb.cellCycleLen + kb.monomerData["degRate"]
			) * monomersView.counts()
		).to('dimensionless').magnitude
```

```nRibosomesNeeded``` can then be compared to the biomass fraction where we know that 15% of the cell is RNA (Dennis and Bremmer), 81% of RNA is rRNA (Neidhardt book), and the weight of rRNA in a single ribosome.

Unfortunatly when using the N-end rule approximatly 2% of the protein degradation rates are ~2 min and the rest ~10 hr. This huge amount of protein turnover caused ```nRibosomesNeeded``` to exceed the theoretical maxiumum set by the biomass function and the data from Dennis, Bremmer, and Neidhardt.

(4) <b>DnaA titration so that initiation occurs at end of cell cycle not beginning</b>
In E. coli replicaiton initiation occurs at the end of the cell cycle (doubling time 60 min) for the next cell. This is the opposite effect from what occured in M. genitalium. DnaA regulates itself - DNA doubling adds more non-specific binding sites which titrates out DnaA so that it stops repressing its own expression - this allows enough DnaA to bind oriC sites.

(5) <b>Relationship between ribosomes and growth - what that means at the molecular level</b> When ribosomes go up what goes down? Does that lead to whole-cell model predictable growth rates? Can we figure this out based on Javier's expression data.

(6) Trade off betwen RNAP and Ribosomes. ppGpp

(7) Show that ribosomes control the growth rate with new FBA implementation with pools (pulling instead of pushin)

(8) Cell start to break down glycogen when it was starving. Dry mass decrease. Don't think any other FBA can predict this.

(9) When elongation rate was increased from ~40 nt/s to ~80 nt/s for rRNA the fraction active of RNAP dropped 5% due to the increase in termination of RNAP on rRNAs. This caused the total rate of growth of major mass fractions to remain relativly invariant. Pretty cool.
