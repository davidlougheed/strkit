# *str*onger

A toolkit for analyzing variation in short(ish) tandem repeats.


## `stronger call`: Genotype caller with bootstrapped confidence intervals

A Gaussian mixture model tandem repeat genotype caller for long read data,
tuned for low-coverage HiFi reads.

### Features:

* Support for re-calling output from `tandem-genotypes`, `RepeatHMM`, and `Straglr`
* Strand resampling / bias correction (for use with the `tandem-genotypes` program)
* 95% confidence intervals on calls via user-configurable bootstrapping

### Notes:

* `--min-allele-reads` will affect the confidence intervals given by the bootstrap process,
  especially in low-coverage loci. This should be set depending on the read technology being used;
  something like a single PacBio HiFi read generally contains higher-quality information than a single
  PacBio CLR read, for example.

### Troubleshooting:

Slow performance can result from running `stronger call` or `stronger re-call` on a system with OpenMP, 
due to a misguided  attempt at multithreading under the hood somewhere in Numpy/Scipy (which doesn't work 
here due to  repeated initializations of the Gaussian mixture model.) To fix this, set the following
environment variable before running:

```bash
export OMP_NUM_THREADS=1
```


## `stronger mi`: Mendelian inheritance analysis

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
