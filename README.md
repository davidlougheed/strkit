# STRkit

A toolkit for analyzing variation in short(ish) tandem repeats.

> **Warning**
> 
> Bootstrapping performance may be hindered on systems with OpenMP without
> additional configuration. See [below](#note-on-openmp-threading) for details
> on how to fix this.


## Installation

STRkit can be installed from PyPI via `pip` with the following command:

```bash
python -m pip install strkit
```


## Commands

### `strkit call`: Genotype caller with bootstrapped confidence intervals

A Gaussian mixture model tandem repeat genotype caller for long read data.
STRkit is tuned specifically for high-fidelity long reads, although other 
long read data should still work.

#### Features:

* Performant, vectorized (thanks to [parasail](https://github.com/jeffdaily/parasail))
  estimates of repeat counts from high-fidelity long reads and a supplied 
  catalog of TR loci and motifs.
* Re-weighting of longer reads, to compensate for their lower likelihood of observation.
  * Whole-genome and targeted genotyping modes to adjust this re-weighting.
* Parallelized for faster computing on clusters and for ad-hoc fast analysis of single samples.
* 95% confidence intervals on calls via a user-configurable optional parametric bootstrapping process.


#### Usage:

```bash
strkit call \
  path/to/read/file.bam \  # [REQUIRED] Indexed read file (BAM/CRAM)
  --ref path/to/reference.fa.gz \  # [REQUIRED] Indexed FASTA-formatted reference genome
  --loci path/to/loci.bed \  # [REQUIRED] TRF-formatted (or 4-col, with motif as last column) list of loci to genotype
  --min-reads 4 \  # Minimum number of supporting reads needed to make a call
  --min-allele-reads 2 \  # Minimum number of supporting reads needed to call a specific allele size 
  --flank-size 70 \  # Size of the flanking region to use on either side of a region to properly anchor reads
```

If you want to output a full call report, you can use the `--json output-file.json` argument to
specify a path to output a more detailed JSON document to. This document contains 99% CIs, peak
labels, and some other information that isn't included in the normal TSV file.

##### Note on OpenMP Threading

Slow performance can result from running `strkit call` or `strkit re-call` on a system with OpenMP, 
due to a misguided  attempt at multithreading under the hood somewhere in Numpy/Scipy (which doesn't work 
here due to  repeated initializations of the Gaussian mixture model.) To fix this, set the following
environment variable before running:

```bash
export OMP_NUM_THREADS=1
```


#### All optional flags:

* `--min-reads ##`: Minimum number of supporting reads needed to make a call. Default: 4
* `--min-allele-reads ##`: Minimum number of supporting reads needed to call a specific allele size. 
  **Default:** 2
* `--flank-size ##`: Size of the flanking region to use on either side of a region to properly anchor reads. 
  **Default:** 70
* `--targeted`: Turn on targeted genotyping mode, which re-weights longer reads differently. Use this option if
  the alignment file contains targeted reads, e.g. from PacBio No-Amp Targeted Sequencing. **Default:** off
* `--num-bootstrap ###`: Now many bootstrap re-samplings to perform. **Default:** 100
* `--sex-chr ??`: Sex chromosome configuration. **Without this, loci in sex chromosomes will not be genotyped.**
  Can be any configuration of Xs and Ys; only count matters. **Default:** *none*
* `--json [path]`: Path to output JSON call data to. JSON call data is more detailed than the `stdout` TSV output.
  **Default:** *none*
* `--no-tsv`: Suppresses TSV output to `stdout`. Without `--json`, no output will be generated, which isn't very 
  helpful. **Default:** TSV output on


### `strkit visualize`: Call visualizer

STRkit bundles a call visualization tool which takes as input a BAM file and
a JSON call file from using the `--json` flag with `strkit call`.

It starts a web server on your local machine; the visualizations can be 
interacted with in a web browser.

To use the tool, run the following command:

```bash
strkit visualize path/to/my-alignment.bam \ 
  --ref hg38 \  # or hg19
  --json path/to/my-calls.json \
  -i 1  # 1-indexed offset in JSON file for locus of interest. Default is 1 if left out.
```

This will output something like the following:

```
 * Serving Flask app 'strkit.viz.server' (lazy loading)
 * Environment: production
   WARNING: This is a development server. Do not use it in a production deployment.
   Use a production WSGI server instead.
 * Debug mode: on
 * Running on http://localhost:5011 (Press CTRL+C to quit)
...
```

You can then go to the URL listed, `http://localhost:5011`, on your local machine
to see the visualization tool:

![Browser Histogram](./docs/images/browser_hist.png)
*STRkit browser histogram, showing an expansion in the HTT gene.*

![igv.js Genome Browser](./docs/images/browser_igv.png)
*The same expansion, shown in the igv.js browser. Note the insertions on
the left-hand side in most reads, and the heterozygous copy number pattern.*

To exit the tool, press `Ctrl-C` in your command line window as mentioned in 
the start-up instructions.


### `strkit re-call`: Genotype re-caller

This command has a similar feature-set as `strkit call`, but is designed to
be used with the output of other long-read STR genotyping methods to refine
the genotype estimates when calling from HiFi reads.\

#### Features:

* Support for re-calling output from `tandem-genotypes`, `RepeatHMM`, and `Straglr`
* Strand resampling / bias correction (for use with the `tandem-genotypes` program)
* 95% confidence intervals on calls via user-configurable bootstrapping

#### Notes:

* `--min-allele-reads` will affect the confidence intervals given by the bootstrap process,
  especially in low-coverage loci. This should be set depending on the read technology being used;
  something like a single PacBio HiFi read generally contains higher-quality information than a single
  PacBio CLR read, for example.



### `strkit mi`: Mendelian inheritance analysis

This tool is currently in development and in a very unfinished state. However, the following features
will be in the final release:

* Mendelian inheritance % (MI) calculations for many common TR genotyping tools for both long/short reads
* Confidence-interval MI calculations for the genotyping tools which report CIs
* Reports of loci (potentially of interest) which do not respect MI



## Copyright and License

Copyright (C) 2021-2022  David Lougheed & McGill University

Portions (`viz`) copyright (C) 2021-2022  David Lougheed

Portions of `viz/templates/browser.html` copyright (C) 2021-2022  Observable, Inc.
Used under the terms of the ISC license.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
