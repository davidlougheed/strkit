# Advanced caller usage and configuration


## All optional flags

* `--sample-id example_sample`: Set a sample ID, or override the alignment file sample ID. This will be included in JSON 
  output, but not TSV output.
* `--min-reads ##`: Minimum number of supporting reads needed to make a call. **Default:** 4
* `--min-allele-reads ##`: Minimum number of supporting reads needed to call a specific allele size. 
  **Default:** 2
* `--max-reads ##`: Maximum number of supporting reads before a locus is skipped. **Default:** 250
* `--min-avg-phred ##`: Minimum average PHRED score for relevant bases (flanking region + tandem repeat).
  Read segments with average PHRED scores below this (common with a threshold of ~13 and ONT Ultra Long reads, 
  for example) will be skipped. **Default:** 13
* `--flank-size ##`: Size of the flanking region to use on either side of a region to properly anchor reads. 
  **Default:** 70
* `--realign` or `-a`: Whether to perform local re-alignment to attempt recovery of soft-clipped reads. Some aligners
  may soft-clip around large insertions, e.g. with an expansion (I've noticed this with *pbmm2*/*minimap2*). 
  Currently recommended **for HiFi only**, since this step aggressively filters out realignments with many mismatches 
  or small indels. Enabling this slows down calling, so it may not be suitable for a very large catalog of STRs.
* `--hq`: Whether to treat provided reads as "high quality", i.e., fairly close to the actual true sequence. Used when 
  detecting expansions, to skip a smoothing filter that may ignore disparate, rare expansion-like read counts.
  Use for CCS reads or similar data (e.g., accurate nanopore sequences) ONLY! **Default:** off
* `--use-hp`: Whether to incorporate `HP` tags from a haplotagged alignment file. This should speed up runtime and 
  will potentially improve calling results. **This flag is experimental, and has not been tested extensively.**
* `--incorporate-snvs [path]` or `--snv [path]`: A path to a VCF with SNVs to incorporate into the calling process and 
  final output. This file is just used as an SNV loci catalog; STRkit itself will perform the SNV calling. Empirically 
  improves calling quality a small amount, speeds up runtime, and gives nearby SNV calls for downstream analysis.
  You can find dbSNP VCFs at
  [`https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/`](https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/).
  The file for GRCh38 is called `00-common_all.vcf.gz` as of time of writing.
* `--targeted` or `-t`: Turn on targeted genotyping mode, which re-weights longer reads differently. Use this option if
  the alignment file contains targeted reads that do not reflect normal mapping patterns, e.g. from PacBio No-Amp 
  Targeted Sequencing. **Default:** off
* `--fractional` or `-f`: Turn on fractional genotyping mode, which allows for partial copy numbers in the reference and 
  in allele calls. *Experimental!* **Default:** off
* `--respect-ref` or `-e`: Turn off reference TR region 'coordinate extension' from what is specified in the catalog.
  TR boundaries can be blurry, so by default we give STRkit an opportunity to extend the provided region to improve
  mapped indel capturing and to be consistent with the approach we use to count repeat copies in non-reference samples.
  Turning this off should give results closer to other STR callers, at the cost of potentially missing variation.
* `--count-kmers` or `-k`: Turn on motif-sized k-mer counting at the allele level, with `-k peak`, or at the read 
  level, with `-k read`, or both with `-k both`. If the flag is provided with no value, it will default to `peak.`
  Note that k-mer counts will only be reported if a `--json` path is specified. This feature can be used to detect
  motif composition differences between alleles or samples. **Default:** `none`
* `--consensus` or `-c`: Turn on consensus calculation for alleles. This adds runtime, but gives a better idea of STR 
  structure and is useful for comparing alleles beyond copy number. If `--vcf` is set, this option is forced on. 
  **Default:** off
* `--num-bootstrap ###` or `-b`: Now many bootstrap re-samplings to perform. **Default:** 100
* `--sex-chr ??` or `-x`: Sex chromosome configuration. **Without this, loci in sex chromosomes will not be genotyped.**
  Can be any configuration of Xs and Ys; only count matters. **Default:** *none*
* `--json [path]` or `-j`: Path to output JSON call data to. JSON call data is more detailed than the `stdout` TSV 
  output. If the value passed is `stdout`, the JSON data will be written to `stdout` instead of a file. 
  **Default:** *none*
* `--indent-json` or `-i`: If passed alongside `--json [x]`, the JSON output will be indented to be more human-readable
  but less compact. **Default:** off
* `--vcf [path]`: Path to output VCF-formatted call data to. Setting this option forces the `--consensus` option as 
  well in order to output true REF/ALT values, which slows down runtime somewhat. If the value passed is `stdout`, the 
  VCF data will be written to `stdout` instead of a file. **Default:** *none*
* `--no-tsv`: Suppresses TSV output to `stdout`. Without `--json`, no output will be generated, which isn't very 
  helpful. **Default:** TSV output on
* `--seed`: Seed the random number generator used for all random sampling, Gaussian mixture modeling, etc. 
  Useful for replicability.


## Usage on HPC machines

We have tested STRkit on three different clusters associated with the 
Digital Research Alliance of Canada (formerly Compute Canada). 

Usage is pretty straightforward; for our use cases we set up a Python virtual environment
with the `strkit` package installed, and ran a SLURM batch job which looks something like:

```bash
#!/bin/bash
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=1-00
#SBATCH --account=rrg-xxxxx


module load python/3.9

cd /home/xxxxx || exit
source env/bin/activate

export OMP_NUM_THREADS=1  # Legacy, should be automatic now but drastically improved performance
strkit call \
  --loci /path/to/catalog \
  --ref /path/to/ref.fa.gz \
  --processes 10 \
  --seed 342 \
  path/to/sample.bam

deactivate

```
