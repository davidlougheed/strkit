# Converting TRF DAT files to BED 

**Version:** v0.25.0+

[Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html) ([Benson 1999](https://pubmed.ncbi.nlm.nih.gov/9862982/)) 
is a widely-used tool for discovering tandem repeat regions in genomes. However, its output is in the form of a `.dat` 
file, which is not always convenient for use in other programs (e.g., TR callers, web browsers)

STRkit provides a built-in subcommand, `strkit convert`, which can (among other things) convert TRF files from the TRF
`.dat` format to either a standard four-column BED file, or a UCSC-style TRF BED file.

To convert a TRF `.dat` file to a sorted four-column BED using STRkit, run the following command:

```bash
strkit convert --in-format trf --out-format bed4 --sort trf.dat > repeats.bed
```

To convert a TRF `.dat` file to a sorted UCSC-style TRF BED file, run the following command instead:

```bash
strkit convert --in-format trf --out-format trf --sort trf.dat > trf.bed
```
