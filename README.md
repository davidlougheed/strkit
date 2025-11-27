# STRkit - short tandem repeat genotyping with long reads

[![PyPI version](https://badge.fury.io/py/strkit.svg)](https://badge.fury.io/py/strkit)
[![BioRxiv DOI](https://img.shields.io/badge/bioRxiv-10.1101/2025.03.25.645269-B31B1B.svg)](https://doi.org/10.1101/2025.03.25.645269)
[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12689906.svg)](https://doi.org/10.5281/zenodo.12689906)

STRkit is a short tandem repeat (STR) genotyping and analysis toolkit for long read sequencing data, especially 
PacBio HiFi data. The STRkit software package is written in Python and is available in the PyPI package registry or as
a Docker container.

If you use STRkit in published work, please cite our preprint:

> [STRkit: precise, read-level genotyping of short tandem repeats using long reads and single-nucleotide variation.](https://doi.org/10.1101/2025.03.25.645269)
> David R Lougheed, Tomi Pastinen, Guillaume Bourque. *BioRxiv&nbsp;preprint*.
> DOI:&nbsp;[10.1101/2025.03.25.645269](https://doi.org/10.1101/2025.03.25.645269)

<img src="./docs/images/strkit_logo_small.png" alt="" width="500" height="324" />

**Officially-supported platforms:** macOS (aarch64), Linux (amd64/aarch64)

> **Note:** Currently, it is unlikely that STRkit will work on "plain" Windows due to downstream dependencies. 
> To run STRkit on Windows, please use [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/) or
> use our Docker image (see [Installation as a Docker container](#as-a-docker-container)).


## Table of Contents

* [Installation](#installation)
  * [Via PyPI](#via-pypi)
  * [As a Docker container](#as-a-docker-container)
* [Commands](#commands)
  * [`strkit call`: Genotype caller with bootstrapped confidence intervals](#strkit-call-genotype-caller-with-bootstrapped-confidence-intervals)
    * [Features](#features)
    * [Usage](#usage)
    * [Further documentation on the STRkit caller, including output format](#further-documentation-on-the-strkit-caller-including-output-format)
  * [`strkit visualize`: Call visualizer](#strkit-visualize-call-visualizer)
  * [`strkit mi`: Mendelian inheritance analysis](#strkit-mi-mendelian-inheritance-analysis)
    * [Usage](#usage-1)
    * [Further documentation](#further-documentation)
  * [`strkit convert`: STR catalog conversion](#strkit-convert-str-catalog-conversion)
    * [Usage](#usage-2)
* [Copyright and License](#copyright-and-license)
  * [Notice](#notice)
  * [Exceptions](#exceptions)


## Installation

### Via PyPI

STRkit requires Python 3.10+ and can be installed from PyPI via `pip` 
with the following command:

```bash
python -m pip install strkit
```

You may need to install the [Rust toolchain](https://www.rust-lang.org/tools/install) (version 1.85.0 or newer)
and a C compiler (e.g., `gcc`, `clang`), as well as `cmake`, to compile the `strkit_rust_ext` wheel, 
although prebuilt wheels for this module are available for some platforms. Compiling the wheel may take quite
a long time (in the tens of minutes).

On Digital Research Alliance of Canada/Compute Canada clusters, this involves loading a few modules:

```bash
module load rust/1.85.0 clang/18.1.8 python/3.11 scipy-stack/2025a parasail/2.6.2
python -m pip install strkit
```

STRkit should then be available in your Python environment as a command-line tool:

```bash
strkit --help
```

### As a Docker container

STRkit is also available as a [Docker container](https://github.com/davidlougheed/strkit/pkgs/container/strkit), stored 
in the GitHub Container Registry.

It can be pulled using the following command:

```bash
docker pull ghcr.io/davidlougheed/strkit:latest
```

Then, STRkit commands can be run mostly as normal using the Docker image:

```bash
docker run -it ghcr.io/davidlougheed/strkit --help
```

For example, to run the `strkit call` subcommand in Docker, including writing file output (via binding input files and 
mounting a read/writeable folder for outputs):

```bash
docker run \
  -v inputs:/inputs:ro \
  -v out:/out \
  -it ghcr.io/davidlougheed/strkit call \
  /inputs/file.bam \
  --hq --loci /inputs/loci.bed --ref /inputs/ref.fa --incorporate-snvs /inputs/dbsnp.vcf.gz \
  --vcf /out/calls.vcf --no-tsv
```

Inside the `inputs` directory, STRkit would (given these parameters) expect the following files to be present:

* `file.bam`: A readset for the sample to be called.
* `file.bam.bai`: A BAI index for the readset.
* `ref.fa`: The reference genome for the reads in FASTA format.
* `ref.fa.fai`: A FAIDX-format index for the reference genome.
* `dbsnp.vcf.gz`: A VCF of SNVs to use as a catalogue for calling SNVs to aid in the STR genotyping process.
* `dbsnp.vcf.gz.tbi`: A Tabix-format index for the dbSNP VCF.
* `loci.bed`: A tab-delimited BED-format catalogue of loci (see the 
  [Caller catalog format & choosing a catalog](./docs/caller_catalog.md) document for more information).

For more information on call parameters, see the below section on `strkit call` and the 
[Advanced caller usage and configuration](./docs/caller_usage.md) document.


## Commands

### `strkit call`: Genotype caller with bootstrapped confidence intervals

A Gaussian mixture model tandem repeat genotype caller for long read data.
STRkit is tuned specifically for high-fidelity long reads, although other 
long read data should still work.

![Calling approach flow chart](./docs/images/call_method_flow.png)

#### Features:

* Performant, vectorized (thanks to [parasail](https://github.com/jeffdaily/parasail))
  estimates of repeat counts from high-fidelity long reads and a supplied 
  catalog of TR loci and motifs.
* Re-weighting of longer reads, to compensate for their lower likelihood of observation.
  * Whole-genome and targeted genotyping modes to adjust this re-weighting.
* Incorporation of single-nucleotide variation (SNVs) for better and faster calling plus 
  additional downstream analysis possibilities.
  * Recommended for **HiFi data and ONT R10 data only**. In my testing, this worsens runtime and call quality for 
    ONT ultra-long-read data, but speeds up the tool and improves call quality for HiFi/ONT R10 data. 
  * **Important note:** This functionality is best for whole-genome STR surveying. If hunting for rare / low-coverage
    expansions, or using targeted sequencing data, it may be best to **NOT USE** an SNV catalog to avoid discording 
    reads with short flanking regions.
* Parallelized for faster computing on clusters and for ad-hoc fast analysis of single samples.
* 95% confidence intervals on calls via a user-configurable optional parametric bootstrapping process.


#### Usage:

See all parameters and example usage with a Slurm cluster: 
[Advanced caller usage and configuration](./docs/caller_usage.md)

##### EXAMPLE USAGE

```bash
# For the dbSNP VCF used below for SNV incorporation, see https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/
# (00-common_all.vcf.gz)
#
# "Accurate reads" here means HiFi / ONT R10 duplex reads, but in practice may also include ONT R10 simplex reads.

strkit call \
  path/to/read/file.bam \  # [REQUIRED] One indexed read file (BAM/CRAM)
  --hq \  # If using accurate reads, enable this to get better genotyping & more robust expansion detection
  --realign \  # If using accurate reads, enable this to enable local realignment / read recovery. Good for detecting expansions, but slows down calling.
  --ref path/to/reference.fa.gz \  # [REQUIRED] Indexed FASTA-formatted reference genome
  --loci path/to/loci.bed \  # [REQUIRED] TRF-formatted (or 4-col, with motif as last column) sorted list of loci to genotype
  --incorporate-snvs path/to/dbsnp/00-common_all.vcf.gz \   # If you want, specify a SNV catalogue to help phase STRs & speed up calling
  --vcf my-calls.vcf \  # Calculate consensus sequences for alleles and output a .vcf (or .vcf.gz) with call data
  --seed 183 \  # Fixed random number generator seed for replicability
  --processes 10 \  # Number of parallel processes to use; DEFAULT: 1
  --no-tsv  # If VCF output is enabled as above, we don't need TSV genotype output to stdout (which is the default)
```

##### REGARDING ALIGNMENTS

Ideally, you should be using a read file aligned with parameters tuned for tandem repeats. 
PacBio provides a 
[recommended workflow](https://github.com/PacificBiosciences/apps-scripts/tree/master/RepeatAnalysisTools)
for CCS alignment in this scenario. However, regular aligned readsets are fine and have been tested
extensively.

If you're using accurate long reads (e.g., HiFi, ONT R10 duplex) as input, **use the `--hq` and 
`--realign` options** to get better genotype calculation and a greater proportion of reads 
incorporated into the computed genotypes, respectively. These should not add much performance 
overhead. *In practice, these options may also aid calling with slightly-less-accurate reads.*

If you want to **incorporate haplotagging from an alignment file (`HP` tags)** into the 
process, which should speed up runtime and potentially improve calling results, you must pass 
the `--use-hp` flag.

##### REGARDING SNV INCORPORATION

If you want to **incorporate SNV calling** into the process, which speeds up runtime and gives
marginally better calling results, you must provide an indexed, `bgzip`-compressed SNV catalog 
VCF which matches your reference genome. You can find dbSNP VCFs at
[`https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/`](https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/).
The file for GRCh38 is called `00-common_all.vcf.gz` as of time of writing.
**Note that this does not need to be an SNV call file for your sample, specifically**; just one 
which has positions, reference/alternate alleles, and the `ID` field populated.

##### REGARDING OUTPUT

If you want to output a full call report, you can use the `--json output-file.json` argument to
specify a path to output a more detailed JSON document to. This document contains 99% CIs, peak
labels, and some other information that isn't included in the normal TSV file. If you want this
file to be indented and human-readable, use the `--indent-json` flag in addition to `--json ...`.

If you want to output a VCF file (STRs and SNVs if called; currently not phased), use the
`--vcf ...` argument. If you pass `--vcf stdout`, the VCF will be written to `stdout` instead of a 
file.

For more information, see also documentation on the [Output formats](./docs/output_formats.md).

###### Note on output locus coordinates

STRkit can adjust reference coordinates of STR loci to include regions it thinks should be part of
the definition. Thus, output coordinates may not match the original BED catalog definition. This is
done for two purposes:
* To include slightly mismapped indels that may lie outside the catalog-defined region
* To try to achieve a more consistent absolute copy number between the reference region and the
  sample data

Using the `--respect-ref` flag (see [Caller usage](./docs/caller_usage.md)), users may suppress 
this behaviour.

##### REGARDING REFERENCE GENOMES

The reference genome provided must be BGZipped and indexed using `samtools faidx`:

```bash
# Starting from a .fa:
bgzip my-reference.fa  # Replaces .fa with a .fa.gz file
samtools faidx my-reference.fa.gz  # Generates a .fai index file
```

##### OTHER PARAMETERS

See the '[Caller catalog format & choosing a catalog](./docs/caller_catalog.md)' page for more on
how to format a locus catalog or choose from existing available catalogs.


#### Further documentation on the STRkit caller, including output format:

  * [Advanced caller usage and configuration](./docs/caller_usage.md)
  * [Caller catalog format & choosing a catalog](./docs/caller_catalog.md)
  * [Output formats](./docs/output_formats.md)


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


#### Further documentation on the STRkit web visualization tool:

* [Using `strkit visualize` remotely](./docs/remote_viz.md)



### `strkit mi`: Mendelian inheritance analysis

Using trio data, candidate de novo STR mutations (or genotyping errors/dropout rates) can be discovered 
by looking at inheritance patterns. This tool provides a few different ways to do this, via:

* Mendelian inheritance % (MI) calculations for many common TR genotyping tools for both long/short reads, 
  including support for genotyping methods which report confidence intervals.
* Reports of loci (potentially of interest) which do not respect MI.

#### Usage

For a basic JSON report on Mendelian inheritance with a trio of STRkit VCFs (compressed and indexed with BGZip), use 
something like the following command:

```bash
# In addition to summary figures on Mendelian inheritance, this tool outputs loci which do not respect MI, which may be 
# useful as candidate de novo mutations. The --mismatch-out-mi flag controls which form of MI metric is used for 
# deciding which loci to output. Options for this flag are:
#   strict (strict copy number MI),
#   pm1 (copy number MI ± 1 repeat unit),
#   ci_95 (copy number 95% confidence interval),
#   ci_99 (copy number 99% confidence interval),
#   seq ([allele] sequence MI),
#   sl ([allele] sequence length MI),
#   sl_pm1 ([allele] sequence length MI ± 1 base pair)
strkit mi \
  --caller strkit-vcf \
  --json mi-report.json \
  --mismatch-out-mi seq \
  child-calls.vcf.gz \
  mother-calls.vcf.gz \
  father-calls.vcf.gz
# This will also output a TSV report to stdout. If this is not desired, use --no-tsv to suppress TSV output.
```

For other options and what they do, run `strkit mi` (with no other arguments) or `strkit mi --help`.

#### Further documentation

**For more information on what kind of analyses can be done with this data**, see the
[Trio analyses with STRkit](./docs/trio_analyses.md) page.


### `strkit convert`: STR catalog conversion

STRkit takes as input a four-or-more-column BED file, structured like:

```
contig  start end [0 or more extraneous columns] motif
```

Any extraneous columns are removed, (internally) leaving a four-column STR locus representation. 
Some other tools, e.g., [Straglr](https://github.com/bcgsc/straglr), also take a four-column STR
BED as locus catalog input. However, other formats representing a catalog of STRs exist:

* [Tandem Repeats Finder](https://github.com/Benson-Genomics-Lab/TRF) outputs a TSV/BED with a lot 
  of information. This can be used as-is with STRkit, but it's safer for other tools to convert to
  a four-column BED format.
* [TRGT uses a custom repeat definition format](https://github.com/PacificBiosciences/trgt/blob/main/docs/repeat_files.md),
  which can specify more advanced STR structures.

#### Usage

The `strkit convert` sub-command requires an input format (`trf` or `trgt`), an output format 
(many, see `strkit convert --help`), and an input file. Output is written to `stdout`.

*Note:* Not all input/output format pairs have available converter functions; an error will be 
printed to `stderr` if one does not exist.

For example, to convert from a TRF BED to a TRGT repeat definition BED file:

```bash
strkit convert --in-format trf --out-format trgt in_file.trf.bed > out_file.bed
```

To attempt a conversion from a TRGT repeat definition file to a STRkit/four-column motif BED:

```bash
strkit convert --in-format trgt --out-format strkit in_file.trgt.bed > out_file.bed
```

Note that TRGT can represent STRs with complex structure that STRkit cannot, so some of these loci
may not be converted (these will be logged to `stderr`).


## Copyright and License

* 2021-2023: &copy; David Lougheed (DL) and McGill University 2021-2023 (versions up to and including `0.8.0a1`), 
  created during graduate research by DL.
* 2023+: (versions beyond `0.8.0a1`):
  * Portions &copy; DL and McGill University 2021-2023
  * Portions &copy; McGill University 2024-2025
  * Portions &copy; DL 2024-2025


### Notice

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

### Exceptions

**Some exclusions to this license apply; specifically portions of 
[`strkit/viz/templates/browser.html`](strkit/viz/templates/browser.html) and
the STRkit logo files ([./docs/images/strkit_logo_small.png](./docs/images/strkit_logo_small.png)
and [./strkit/viz/static/logo.png](./strkit/viz/static/logo.png).)**

The STRkit logo is &copy; David Lougheed 2022, and was designed by Evelyn Lougheed. It is not licensed
under the terms of the GPL 3.0; it is instead licensed under the terms of the 
[CC BY-ND 4.0](https://creativecommons.org/licenses/by-nd/4.0/).

Portions of `viz/templates/browser.html` copyright (C) 2021-2022  Observable, Inc.
Used under the terms of the ISC license.
