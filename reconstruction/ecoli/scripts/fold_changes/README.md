Explore data source of foldChanges.tsv in reconstruction and allow replacement of suspect values.

File sources:
- kb.tsv: output/data EcoliFC/KB.txt (from EcoliFoldChanges repo - 8580aae)
- shifts.tsv: data/EcoMAC/experimentalShift.tsv (from EcoliFoldChanges repo - 8580aae)
- gene_names.tsv: data/EcoMAC/geneName.tsv (from EcoliFoldChanges repo - 8580aae)
- fold_changes.tsv: reconstruction/ecoli/flat/foldChanges.tsv (from wcEcoli repo - 143cf74856)
- updated_fold_change.tsv: generated with `./calculate_fc.py --replace`

Running:
- `./calculate_fc.py`: attempts to use best information available (regulation direction discrepancies)
- `./calculate_fc.py -v`: verbosely print mismatched values and a summary of each TF
- `./calculate_fc.py -m`: attempts to match expected processing of source data (relies on unknown method of determining regulation direction after processing and compares based on assuming only positive direction)
- `./calculate_fc.py --replace`: replaces inconsistent fold changes with new means

Notes:
Original data appears to ignore direction of regulation in data analysis.
Mean and std are calculated based on abs of fold change and direction is
probably set afterwards based on curated regulation.  This can lead to issues
with ambiguous regulation where annotation point to positive and negative
regulation (eg ArcA -> sdhCDAB is positive even though fold change data points
to a negative regulation).

\>90% of data appears to be consistent (82/938 fold changes do not match with -m option and ignoring fold changes missing in wcm).
