# STRkit output formats

STRkit can output three different file formats, depending on the set of arguments used:

* [TSV](#tsv-standard-output): by default, printed to `stdout` when STRkit is run. Good as an overview, but less 
  informative/interoperable than other formats.
* [JSON](#json-report): a JSON report, containing the maximum amount of information possible. These files can be quite 
  large, especially if formatted to be human-readable and indented with the `--indent-json` flag.
* [VCF](#vcf): a [VCF 4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file, with STR and SNV genotypes, including
  consensus STR sequences.

**Note:** In general, the JSON format contains the most information about how STRkit was run, and each locus' called 
genotype.


## TSV (standard output)

A tab-separated text file with the following columns:

* Chromosome
* Starting position (matching input BED file; real coordinates of region may be different if 
  `--respect-ref` is not used)
* Ending position (matching input BED file; real coordinates of region may be different if 
  `--respect-ref` is not used)
* Motif sequence (matching input BED file)
* Reference copy number
* Comma-delimited list of copy numbers for all reads successfully extracted for this locus.
* Copy number call, `|`-delimited (one call per allele)
* 95% confidence intervals for copy number calls, `|`-delimited (one `X-Y` 95% CI per allele)
* Calling approach used by STRkit: one of:
  * `dist` - clustering based on a copy number distance metric
  * `snv+dist` - clustering based on a copy number + nearby SNV genotype difference distance metric
  * `snv` - clustering solely based on nearby SNV genotypes

Here is an example line:

```
chr4	5975495	5975530	TTTTG	7	6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8	6|7	6-6|7-7	snv
```

Note that quite a bit of information is missing from the TSV


## JSON report

TODO


## VCF

TODO
