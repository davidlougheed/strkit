# *str*kit

A toolkit for analyzing variation in short(ish) tandem repeats.


## `strkit call`: Genotype caller with bootstrapped confidence intervals

A Gaussian mixture model tandem repeat genotype caller for long read data,
tuned for high-fidelity long reads.

### Features:

* Performant, vectorized (via [parasail](https://github.com/jeffdaily/parasail))
  estimates of repeat counts from high-fidelity long reads and a supplied 
  catalog of TR loci.
* Parallelized for faster computing on clusters.
* 95% confidence intervals on calls via user-configurable bootstrapping.


### Usage:

```bash
strkit call \
  path/to/read/file.bam \  # [REQUIRED] Indexed read file (BAM/CRAM)
  --ref path/to/reference.fa.gz \  # [REQUIRED] Indexed FASTA-formatted reference genome
  --loci path/to/loci.bed \  # [REQUIRED] TRF-formatted (or 4-col, with motif as last column) list of loci to genotype
  --min-reads 4 \  # Minimum number of supporting reads needed to make a call
  --min-allele-reads 2 \  # Minimum number of supporting reads needed to call a specific allele size 
  --flank-size 70 \  # Size of the flanking region to use on either side of a region to properly anchor reads
```

Slow performance can result from running `strkit call` or `strkit re-call` on a system with OpenMP, 
due to a misguided  attempt at multithreading under the hood somewhere in Numpy/Scipy (which doesn't work 
here due to  repeated initializations of the Gaussian mixture model.) To fix this, set the following
environment variable before running:

```bash
export OMP_NUM_THREADS=1
```



## `strkit re-call`: Genotype re-caller

This command has the same feature-set as `strkit call`, but is designed to
be used with the output of other long-read STR genotyping methods to refine
the genotype estimates when calling from HiFi reads.\

### Features:

* Support for re-calling output from `tandem-genotypes`, `RepeatHMM`, and `Straglr`
* Strand resampling / bias correction (for use with the `tandem-genotypes` program)
* 95% confidence intervals on calls via user-configurable bootstrapping

### Notes:

* `--min-allele-reads` will affect the confidence intervals given by the bootstrap process,
  especially in low-coverage loci. This should be set depending on the read technology being used;
  something like a single PacBio HiFi read generally contains higher-quality information than a single
  PacBio CLR read, for example.



## `strkit mi`: Mendelian inheritance analysis

This tool is currently in development and in a very unfinished state. However, the following features
will be in the final release:

* Mendelian inheritance % (MI) calculations for many common TR genotyping tools for both long/short reads
* Confidence-interval MI calculations for the genotyping tools which report CIs
* Reports of loci (potentially of interest) which do not respect MI



## Copyright and License

Copyright (C) 2021-2022  David Lougheed & McGill University

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
