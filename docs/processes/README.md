## Processes

The Processes of the *E. coli* model span several major areas of cellular physiology. We have clustered the Processes into the following groups, which we present in the order listed: Central Dogma, Metabolism, and Balanced Growth. We modeled Processes using the most appropriate mathematics for their individual network topology and degree of experimental characterization. Each process is a computational representation of chemical reactions or transformations grouped by a physiological function. The actual division of reactions across processes is a modeling decision made during model construction, and the number of Processes does not reflect their complexity or scope. The inputs and outputs of each Process are the counts of metabolites or macromolecules and the catalytic capacity or configuration of the enzymes that catalyze the reactions in each Process. This section details the model implementation, computational algorithm, associated data, and relevant code for each Process.

### Process Index

* [Cell Division](cell_division.pdf)
* [Chromosome Replication](chromosome_replication.pdf)
* [Complexation](complexation.pdf)
* [Metabolism](metabolism.pdf)
* [Protein Degradation](protein_degradation.pdf)
* [RNA Degradation](rna_degradation.pdf)
* [Transcription](transcription.pdf)
* [Transcripton Regulation](transcription_regulation.pdf)
* [Translation](translation.pdf)

### Regeneration of pdf files from tex files

Assuming all the dependencies have been installed, the pdf files can be compiled from the tex files by changing to the `src` directory and running the following commands

      pdflatex PROCESS_NAME.tex
      bibtex PROCESS_NAME.aux
      pdflatex PROCESS_NAME.tex
      pdflatex PROCESS_NAME.tex
      mv PROCESS_NAME.pdf ../

where `PROCESS_NAME` is replaced by the name of the process for which the `tex` file has been edited.
