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

Note that quite a bit of information is missing from the TSV, including per-sample copy numbers, read identities, 
SNV calls, and STR consensus sequences.


## JSON report

Example report format:

```javascript
{
  "sample_id": "HG002",
  "caller": {
    "name": "strkit",
    "version": "0.15.0"
  },
  "parameters": {
    "read_files": "HG002.SequelII.ccs.phased.40x.chr4.bam",
    "reference_file": "/Users/davidlougheed/git/gt-poc/hg38.analysisSet.fa.gz",
    "min_reads": 4,
    "min_allele_reads": 2,
    "min_avg_phred": 13,
    "num_bootstrap": 100,
    "flank_size": 70,
    "sample_id": "HG002",
    "realign": true,
    "hq": true,
    "snv_vcf": "00-common_all.vcf.gz",
    "snv_min_base_qual": 20,
    "targeted": false,
    "respect_ref": false,
    "count_kmers": "none",
    "consensus": true,
    "log_level": 10,
    "seed": 1234,
    "processes": 1
  },
  "runtime": 8.628772,
  "contigs": [
    "chr4"
  ],
  "results": [
    {
      "locus_index": 1,
      "contig": "chr4",
      "start": 96617,
      "end": 96648,
      "start_adj": 96617,
      "end_adj": 96648,
      "motif": "AC",
      "ref_cn": 16,
      "ref_start_anchor": "t",
      "ref_seq": "acacacacacacacacacacacacacacaca",
      "reads": {
        "m64011_190901_095311/50792740/ccs": {
          "s": "-",
          "cn": 15,
          "w": 1.0217145751733625,
          "snvu": ["G"],
          "p": 0
        },
        // ...
        "m64012_190921_234837/4523939/ccs": {
          "s": "+",
          "cn": 15,
          "w": 1.0217145751733625,
          "snvu": ["A"],
          "p": 1
        },
        // ...
      },
      "snvs": [
        {
          "id": "rs73213545",
          "ref": "G",
          "pos": 94593,
          "call": ["G", "A"],
          "rcs": [20, 23]
        }
      ],
      "assign_method": "snv+dist",
      "call": [15, 15],
      "call_95_cis": [
        [15, 15],
        [15, 15]
      ],
      "call_99_cis": [
        [15, 15],
        [15, 15]
      ],
      "peaks": {
        "means": [15, 15],
        "weights": [0.5, 0.5],
        "stdevs": [0.31622776601683794, 0.3585309239667531],
        "modal_n": 2,
        "n_reads": [20, 23],
        "seqs": [
          ["ACACACACACACACACACACACACACACA", "poa"],
          ["ACACACACACACACACACACACACACACA", "poa"]
        ]
      },
      "read_peaks_called": true,
      "time": 0.1274
    },
    // ...
  ]
}
```


## VCF

VCF format fields (i.e., for each variant sample entry):

* `AD`: Read depth for each allele
* `DP`: Total read depth
* `DPS`: Total read depth; only supporting reads (for calls with incorporated SNVs mainly; STR calls only)
* `GT`: Genotype
* `MC`: Motif copy number for each allele (STR calls only)
* `MCCI`: Motif copy number 95% confidence intervals for each allele (STR calls only)
* `PS`: Phase set
* `PM`: Peak-calling method (`dist`/`single`/`snv+dist`/`snv`/`hp`; STR calls only)

VCF info. fields (i.e., for each STR variant record; not present for SNV records):

* `VT`: Variant record type (`str` or `snv`)
* `MOTIF`: Motif sequence
* `REFMC`: Motif copy number in the reference genome
